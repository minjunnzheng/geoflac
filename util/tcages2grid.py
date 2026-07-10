#!/usr/bin/env python3

'''Average marker thermochron ages onto grid elements and save them in grid VTS files.

Each element value is the arithmetic mean over all markers inside the element,
with open-system markers counted as age 0 and unreset markers counted as
age = model_time + 100 Myr (the unreset_time sentinel of the Fortran
geoflac-thermochron implementation), reproducing its smooth chronage grids.
'''

from __future__ import print_function
import sys
import os
import glob
import xml.etree.ElementTree as ET
import numpy as np
import io
import flac
import flac2vtk

# Fortran geoflac-thermochron: unreset_time = -100 Myr, so an unreset marker
# contributes age = model_time + 100 Myr to the element average.
unreset_age_offset = 100.0

def read_data_array(el):
    dtype_str = el.get("type", "Float32")
    if dtype_str == "Float32":
        dtype = np.float32
    elif dtype_str == "Int32":
        dtype = np.int32
    else:
        dtype = np.float32 # fallback
        
    fmt = el.get("format", "ascii")
    
    text = el.text
    if text is None:
        text = ""
        
    if fmt == "ascii":
        return np.fromstring(text, sep=' ', dtype=dtype)
    elif fmt == "binary":
        text = text.strip()
        if not text:
            return np.array([], dtype=dtype)
        uncompressed_data = flac.FlacFromVTK._unpack_vtk(None, text)
        return np.frombuffer(uncompressed_data, dtype=dtype).copy()
    else:
        raise ValueError(f"Unknown format: {fmt}")


def make_vts_dataarray_element(data, name, format_type):
    orig_output_in_binary = flac2vtk.output_in_binary
    flac2vtk.output_in_binary = (format_type == "binary")
    
    f_str = io.StringIO()
    flac2vtk.vts_dataarray(f_str, data, name)
    
    flac2vtk.output_in_binary = orig_output_in_binary
    
    xml_fragment = f_str.getvalue()
    return ET.fromstring(xml_fragment)


def interpolate_vtp_to_vts(vtp_path, vts_path):
    print(f"Interpolating {os.path.basename(vtp_path)} -> {os.path.basename(vts_path)}   ", end='\r')
    sys.stdout.flush()
    
    # 1. Parse VTS to get grid coordinates and dimensions
    vts_tree = ET.parse(vts_path)
    vts_root = vts_tree.getroot()
    
    grid_el = vts_root.find(".//StructuredGrid")
    extent_str = grid_el.get("WholeExtent")
    extent_parts = [int(x) for x in extent_str.split()]
    nex = extent_parts[1]
    nez = extent_parts[3]
    nx = nex + 1
    nz = nez + 1
    
    # Check format of VTS data arrays
    first_da = vts_root.find(".//PointData/DataArray")
    if first_da is not None:
        format_type = first_da.get("format", "ascii")
    else:
        format_type = "binary" # fallback
        
    # 2. Parse VTP to get marker coordinates and thermochron ages
    vtp_tree = ET.parse(vtp_path)
    vtp_root = vtp_tree.getroot()
    
    piece_el = vtp_root.find(".//Piece")
    nmarkers = int(piece_el.get("NumberOfPoints", 0))
    
    if nmarkers == 0:
        # If there are no markers, we write all NaNs for any system found in PointData
        point_data_el = vtp_root.find(".//PointData")
        for da in point_data_el.findall("DataArray"):
            name = da.get("Name")
            if name and (name.startswith("age_") or name == "max_temp"):
                elem_age = np.full((nez, nex), np.nan, dtype=np.float32)
                new_da = make_vts_dataarray_element(elem_age, name, format_type)
                
                cell_data_el = vts_root.find(".//CellData")
                found = False
                for cda in cell_data_el.findall("DataArray"):
                    if cda.get("Name") == name:
                        idx = list(cell_data_el).index(cda)
                        cell_data_el[idx] = new_da
                        found = True
                        break
                if not found:
                    cell_data_el.append(new_da)
                    
        # Write VTS back
        xml_str = ET.tostring(vts_root, encoding='utf-8').decode('utf-8')
        with open(vts_path, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write(xml_str)
        return

    # Model time (Myr) of this frame, needed for the unreset sentinel age
    time_el = vtp_root.find(".//FieldData/DataArray[@Name='TIME']")
    if time_el is None:
        raise ValueError("TIME field data not found in VTP!")
    t_frame = float(time_el.text.strip())

    # Extract ntriag
    point_data_el = vtp_root.find(".//PointData")
    ntriag_el = None
    for da in point_data_el.findall("DataArray"):
        if da.get("Name") == "ntriag":
            ntriag_el = da
            break
            
    if ntriag_el is None:
        raise ValueError("ntriag array not found in VTP PointData!")
        
    ntriag_markers = read_data_array(ntriag_el)
    
    # Compute element index for each marker
    ntriag_clean = np.maximum(1, ntriag_markers)
    k_m = (ntriag_clean - 1) % 2 + 1
    j_el = np.clip(((ntriag_clean - k_m) // 2) % (nz - 1), 0, nz - 2)
    i_el = np.clip(((ntriag_clean - k_m) // 2) // (nz - 1), 0, nx - 2)

    # Find all thermochronology systems (starting with age_ or max_temp)
    cell_data_el = vts_root.find(".//CellData")
    for da in point_data_el.findall("DataArray"):
        name = da.get("Name")
        if name and (name.startswith("age_") or name == "max_temp"):
            ages = read_data_array(da)

            # Same scheme as out_thermochron_grid in the Fortran geoflac-thermochron:
            # plain per-element mean over all markers, with the marker sentinels mapped
            # to their Fortran-equivalent numeric ages first.
            vals = ages.astype(np.float64)
            if name != "max_temp":
                vals[ages == -1.0] = 0.0                            # open system: closure time ~ now
                vals[np.isnan(ages)] = t_frame + unreset_age_offset  # unreset

            sum2d = np.zeros((nez, nex))
            cnt2d = np.zeros((nez, nex))
            np.add.at(sum2d, (j_el, i_el), vals)
            np.add.at(cnt2d, (j_el, i_el), 1.0)

            # Elements without markers are NaN (the Fortran code writes 0.0 there,
            # but 0.0 would be indistinguishable from a real zero age).
            elem_age = np.full((nez, nex), np.nan, dtype=np.float32)
            has_marker = cnt2d > 0
            elem_age[has_marker] = (sum2d[has_marker] / cnt2d[has_marker]).astype(np.float32)

            new_da = make_vts_dataarray_element(elem_age, name, format_type)
            
            # Update/Insert into CellData of VTS
            found = False
            for cda in cell_data_el.findall("DataArray"):
                if cda.get("Name") == name:
                    idx = list(cell_data_el).index(cda)
                    cell_data_el[idx] = new_da
                    found = True
                    break
            if not found:
                cell_data_el.append(new_da)
                
    # Write VTS back
    xml_str = ET.tostring(vts_root, encoding='utf-8').decode('utf-8')
    with open(vts_path, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write(xml_str)


def main(args):
    doc = '''usage: tcages2grid.py path [frame_min [frame_max]]
       tcages2grid.py vtp_file vts_file

Average marker thermochron ages onto grid elements and save them in grid VTS files.
'''
    if len(args) < 1:
        print(doc)
        sys.exit(1)
        
    path_or_file = args[0]
    
    if path_or_file.endswith('.vtp'):
        # Mode: Single VTP / VTS pair
        vtp_file = path_or_file
        if len(args) >= 2:
            vts_file = args[1]
        else:
            # Infer VTS file name
            # Replace flacmarker with flac and .vtp with .vts
            base = os.path.basename(vtp_file)
            dirname = os.path.dirname(vtp_file)
            if base.startswith("flacmarker."):
                vts_base = base.replace("flacmarker.", "flac.", 1).replace(".vtp", ".vts", 1)
            else:
                vts_base = base.replace(".vtp", ".vts", 1)
            vts_file = os.path.join(dirname, vts_base)
            
        if not os.path.exists(vtp_file):
            print(f"Error: {vtp_file} not found!")
            sys.exit(1)
        if not os.path.exists(vts_file):
            print(f"Error: {vts_file} not found!")
            sys.exit(1)
            
        interpolate_vtp_to_vts(vtp_file, vts_file)
        print()
        
    else:
        # Mode: Directory with start/end frames
        path = path_or_file
        start = 1
        end = -1
        if len(args) >= 2:
            start = int(args[1])
            if len(args) >= 3:
                end = int(args[2])
                
        # Find matching VTP files
        vtplist = sorted(glob.glob(os.path.join(path, 'flacmarker.*.vtp')))
        if not vtplist:
            print(f"Error: No flacmarker.*.vtp files found in {path}")
            sys.exit(1)
            
        if end == -1:
            end = int(os.path.basename(vtplist[-1])[11:-4])
            
        for i in range(start, end + 1):
            vtp_file = os.path.join(path, f'flacmarker.{i:06d}.vtp')
            vts_file = os.path.join(path, f'flac.{i:06d}.vts')
            if os.path.exists(vtp_file) and os.path.exists(vts_file):
                interpolate_vtp_to_vts(vtp_file, vts_file)
            else:
                print(f"Skipping frame {i}: VTP or VTS file not found.")
        print()

if __name__ == '__main__':
    main(sys.argv[1:])
