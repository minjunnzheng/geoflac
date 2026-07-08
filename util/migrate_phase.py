#!/usr/bin/env python3
"""migrate_phase.py -- deterministic phase-ID migration for geoflac .inp files.

Remaps phase IDs in the nzone_age and inhom blocks (e.g. 13->20, 14->21, ...)
and widens the rheology table so every new ID is a verbatim copy of its source
row; the old rows are left untouched (they become unused). Written for the
`orogen` branch, where orogeny lithologies live in the rule-free phase range
20-24 and change_phase.f90 is kept upstream-original, so a static lithology
that borrows a slot with a built-in conversion rule (16 khydmant, 17 kmetased)
must be migrated before it will run correctly (ARCHITECTURE 3.2/3.4/4.6).

Parser semantics mirror the verified /Users/br/tools/flacdraw/inp.py
(tokens / next_data_line, old & new nzone_age formats, inhom line tracking),
vendored here so the tool is self-contained: python3 only, no dependencies,
runnable locally or on the cluster.

Usage:
  migrate_phase.py in.inp --map "13:20,14:21,16:22,18:23,19:24" -o out.inp
"""
import argparse
import re
import sys
from pathlib import Path

MAXPH = 30  # must match src/params.f90 maxph


# --------------------------------------------------------------------------
# tokenizer (verbatim from flacdraw/inp.py)
# --------------------------------------------------------------------------
def strip_comment(line):
    return line.split(';', 1)[0].strip()


def tokens(line):
    s = strip_comment(line)
    return [t for t in re.split(r'[,\s]+', s) if t] if s else []


def find_section(lines, key, after=0):
    for i in range(after, len(lines)):
        if key in lines[i]:
            return i
    return -1


def next_data(lines, start):
    i = start
    while i < len(lines):
        if tokens(lines[i]):
            return i
        i += 1
    return -1


def trailing_comment(line):
    return ';' + line.split(';', 1)[1] if ';' in line else ''


def leading_ws(line):
    return line[:len(line) - len(line.lstrip())]


# --------------------------------------------------------------------------
# structural parse (indices only -- we rewrite lines, not a data model)
# --------------------------------------------------------------------------
def parse_zones(lines):
    """nzone_age block. Returns (records, index_after_block).
    Each record: {'phase_line': idx, 'fmt': 'old'|'new', 'phases': [int],
                  'tokens': [...] (old only)}."""
    i = find_section(lines, 'nzone_age')
    if i < 0 or not tokens(lines[i]):
        i = find_section(lines, 'Initially Variable thermal age')
        i = next_data(lines, i + 1)
    if i < 0:
        raise ValueError('nzone_age block not found')
    n = int(tokens(lines[i])[0])
    j = next_data(lines, i + 1)
    fmt = 'old' if len(tokens(lines[j])) >= 12 else 'new'
    recs = []
    for _ in range(n):
        if fmt == 'old':
            t = tokens(lines[j])
            recs.append({'phase_line': j, 'fmt': 'old',
                         'phases': [int(float(x)) for x in t[5:10]],
                         'tokens': t})
            j = next_data(lines, j + 1)
        else:
            j2 = next_data(lines, j + 1)           # nph + hc line
            nph = int(tokens(lines[j2])[0])
            j3 = next_data(lines, j2 + 1)          # phase line
            recs.append({'phase_line': j3, 'fmt': 'new',
                         'phases': [int(x) for x in tokens(lines[j3])[:nph]]})
            j = next_data(lines, j3 + 1)
    return recs, j


def parse_rheology(lines, start):
    """Returns (count_idx, nphase, row_indices)."""
    count_idx = next_data(lines, start)
    nphase = int(tokens(lines[count_idx])[0])
    rows, j = [], count_idx
    while len(rows) < nphase:
        j = next_data(lines, j + 1)
        if j < 0:
            raise ValueError('rheology table truncated')
        rows.append(j)
    return count_idx, nphase, rows


def parse_inhom(lines):
    """Returns (count_idx, n, rows[(idx, phase)], irphase)."""
    i = find_section(lines, 'inhom(number')
    if i < 0:
        i = find_section(lines, 'inhom (number')
    irphase = 0
    if i < 0 or not tokens(lines[i]):
        k = find_section(lines, 'irphase')
        irp = next_data(lines, k)
        irphase = int(tokens(lines[irp])[0])
        pf = next_data(lines, irp + 1)          # phasefile / dummy line
        i = next_data(lines, pf + 1)            # inhom count
    count_idx = i
    n = int(tokens(lines[count_idx])[0])
    rows, j = [], count_idx + 1
    while len(rows) < n and j < len(lines):
        t = tokens(lines[j])
        if len(t) >= 7:
            rows.append((j, int(t[4])))
        j += 1
    if len(rows) < n:
        raise ValueError(f'inhom count says {n}, found {len(rows)}')
    return count_idx, n, rows, irphase


# --------------------------------------------------------------------------
# line rewriters (values change, formatting is regenerated but semantics keep)
# --------------------------------------------------------------------------
def rewrite_zone_phase_line(orig, rec, new_ph):
    if rec['fmt'] == 'old':
        t = list(rec['tokens'])
        t[5:10] = [str(p) for p in new_ph]
        return leading_ws(orig) + ', '.join(t) + ' ' + trailing_comment(orig)
    body = leading_ws(orig) + ' '.join(str(p) for p in new_ph)
    tc = trailing_comment(orig)
    return body + ('  ' + tc if tc else '')


def rewrite_inhom_phase_line(orig, new_phase):
    t = tokens(orig)
    t[4] = str(new_phase)
    tc = trailing_comment(orig)
    line = ('  {:>4} {:>4} {:>4} {:>4} {:>5} {:>4} {:>6}'
            .format(*t[:7]))
    return line + ('  ' + tc if tc else '')


# --------------------------------------------------------------------------
# migration
# --------------------------------------------------------------------------
def parse_map(spec):
    mapping = {}
    for pair in spec.split(','):
        pair = pair.strip()
        if not pair:
            continue
        a, b = pair.split(':')
        a, b = int(a), int(b)
        if a in mapping:
            raise SystemExit(f'error: duplicate source id {a} in --map')
        mapping[a] = b
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit('error: --map targets are not unique')
    return mapping


def migrate(text, mapping):
    lines = text.splitlines()
    zrecs, after_zones = parse_zones(lines)
    rcount_idx, nphase, rrows = parse_rheology(lines, after_zones)
    icount_idx, ninhom, irows, irphase = parse_inhom(lines)

    # --- validate the mapping against the actual table -------------------
    if any(k < 1 for k in mapping):
        raise SystemExit('error: source phase ids must be >= 1 (never map -1)')
    for k in mapping:
        if k > nphase:
            raise SystemExit(f'error: source phase {k} > nphase ({nphase}); '
                             'no rheology row to copy')
    for v in mapping.values():
        if v <= nphase:
            raise SystemExit(f'error: target phase {v} <= nphase ({nphase}); '
                             'would overwrite an existing row -- pick a new slot')
        if v > MAXPH:
            raise SystemExit(f'error: target phase {v} exceeds maxph ({MAXPH})')
    if irphase > 0:
        raise SystemExit('error: irphase > 0 (phase distribution read from a '
                         'phasefile); those phases cannot be migrated '
                         'automatically -- handle the phasefile manually')

    inverse = {v: k for k, v in mapping.items()}
    repl = {}          # line index -> new text
    counts = {}        # source id -> number of column cells remapped

    # zones
    for rec in zrecs:
        new_ph = [mapping.get(p, p) for p in rec['phases']]
        for a, b in zip(rec['phases'], new_ph):
            if a != b:
                counts[a] = counts.get(a, 0) + 1
        if new_ph != rec['phases']:
            repl[rec['phase_line']] = rewrite_zone_phase_line(
                lines[rec['phase_line']], rec, new_ph)

    # inhom (phase == -1 is a strain-only seed, never in mapping)
    for idx, ph in irows:
        if ph in mapping:
            counts[ph] = counts.get(ph, 0) + 1
            repl[idx] = rewrite_inhom_phase_line(lines[idx], mapping[ph])

    # rheology: widen to max target, new rows copy their source verbatim
    M = max(mapping.values())
    new_rows = []
    for pid in range(nphase + 1, M + 1):
        src = inverse.get(pid)
        if src is not None:
            new_rows.append(f';({pid}) migrated from ({src})')
            new_rows.append(lines[rrows[src - 1]])
        else:
            new_rows.append(f';({pid}) unused filler (copy of row 1)')
            new_rows.append(lines[rrows[0]])
    repl[rcount_idx] = re.sub(r'^\s*\d+', str(M), lines[rcount_idx], count=1)

    # assemble (in-place replacements + rheology rows spliced after last row)
    last_rheo = rrows[-1]
    out = []
    for idx, ln in enumerate(lines):
        out.append(repl.get(idx, ln))
        if idx == last_rheo:
            out.extend(new_rows)
    return '\n'.join(out) + '\n', counts, nphase, M


def self_verify(text, mapping):
    """Re-parse the migrated text and assert no source id survives in any
    phase column. Returns list of (where, id) violations."""
    lines = text.splitlines()
    zrecs, after_zones = parse_zones(lines)
    _, _, irows, _ = parse_inhom(lines)
    bad = []
    srcs = set(mapping)
    for rec in zrecs:
        for p in rec['phases']:
            if p in srcs:
                bad.append((f'nzone_age line {rec["phase_line"] + 1}', p))
    for idx, ph in irows:
        if ph in srcs:
            bad.append((f'inhom line {idx + 1}', ph))
    return bad


def advisory_checks(text):
    """Best-effort, non-fatal reminders (ARCHITECTURE 2.3 preconditions)."""
    notes = []
    # nofside == 5 is an internal BC that the orogeny setup does not use; in
    # new-format files the BC rows carry no inline keyword, so this is a scan.
    for ln in text.splitlines():
        t = tokens(ln)
        if len(t) >= 5 and 'nofside' in ln.lower() and '5' in t:
            notes.append(f'nofside may be 5: "{ln.strip()}"')
    notes.append('magma extrusion and nofside=5 were NOT auto-verified from a '
                 'new-format .inp (no inline keywords). Confirm magma is OFF '
                 '(ratio_mantle_mzone = 1) and no BC uses nofside=5 -- '
                 'ARCHITECTURE 2.3 preconditions.')
    return notes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inp', help='input .inp file')
    ap.add_argument('--map', required=True,
                    help='comma list old:new, e.g. "13:20,14:21,16:22,18:23,19:24"')
    ap.add_argument('-o', '--out', required=True, help='output .inp file')
    args = ap.parse_args()

    mapping = parse_map(args.map)
    text = Path(args.inp).read_text()
    migrated, counts, nphase, M = migrate(text, mapping)

    bad = self_verify(migrated, mapping)
    if bad:
        for where, pid in bad:
            print(f'  SELF-VERIFY FAIL: source phase {pid} still at {where}',
                  file=sys.stderr)
        raise SystemExit('error: migration did not remove all source ids')

    Path(args.out).write_text(migrated)

    print(f'migrated {args.inp} -> {args.out}')
    print(f'  rheology table widened: nphase {nphase} -> {M}')
    print('  mapping:', ', '.join(f'{k}->{v}' for k, v in sorted(mapping.items())))
    for k in sorted(mapping):
        print(f'    phase {k} -> {mapping[k]}: {counts.get(k, 0)} cells replaced')
    print('  self-verify: OK (no source id remains in any phase column)')
    for note in advisory_checks(text):
        print('  NOTE:', note)


if __name__ == '__main__':
    main()
