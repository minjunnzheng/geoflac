subroutine fl_thermchron
!For Thermochronology
!
!Chase Shyu <iamhemry@gmail.com>
!Jan. 14th, 2013
!Refactored for geoflac-master Jan 2026

use arrays
use marker_data
use params
implicit none

!$ACC routine(temp2marker) seq
character :: screen*200
integer :: kk, i, i_closure_if, nfp, kc
double precision :: t_temp, t_rate, t_closure, delta_temp, sterr
double precision :: t_temp_last, trate_abs, rr1, rr2, cc1, cc2

sterr = 0.d-4
!12 format(1I10,5F15.7)

!if (ireset.eq.1) then
!  iwait = 1
!  return
!end if

if (mod(nloop, 10).ne.0) return
if (i_prestress .eq. 1 ) return

! Update present temperature on marker
!$OMP Parallel do default(none) shared(nmarkers,mark_dead) private(kk)
!$ACC parallel loop async(1)
do kk = 1, nmarkers
  if (mark_dead(kk).eq.0) cycle
  call temp2marker(kk)
end do
!$OMP end parallel do

! do not calculate the closure temperature just after remeshing
if (iwait .eq. 1) then
  iwait = 0
  return
end if

do i = 1, nchron

!$OMP Parallel do default(none) private(kk,t_temp,t_rate,t_temp_last, &
!$OMP                           i_closure_if,t_closure,screen,delta_temp) &
!$OMP                           shared(mark_dead,mark_chron_if,mark_temp,mark_cooling_rate, &
!$OMP                           mark_chron_temp,mark_chron_time, &
!$OMP                           chron_ref,time,dt,i,nloop, &
!$OMP                           sterr,nmarkers,nchron_fpair)
!$ACC parallel loop async(1)
! calculate if closure and assign age
    do kk = 1, nmarkers
      if (mark_dead(kk).eq.0) cycle

      i_closure_if = mark_chron_if(i, kk)

!      if (iupdate_temp_rate .le. 0) cycle
      t_temp = mark_temp(kk)
      t_rate = mark_cooling_rate(kk)
      ! find the closure temperature of t_rate
      if(t_rate .lt. 0.d0)then
        ! inlined closure-temp interpolation (linear scan; device-friendly)
        trate_abs = -1.d0*t_rate
        nfp = nchron_fpair(i)
        if (trate_abs <= chron_ref(i,1,1)) then
          t_closure = chron_ref(i,1,2)
        else if (trate_abs >= chron_ref(i,nfp,1)) then
          t_closure = chron_ref(i,nfp,2)
        else
          t_closure = chron_ref(i,nfp,2)
          do kc = 1, nfp-1
            rr1 = chron_ref(i,kc,1); rr2 = chron_ref(i,kc+1,1)
            if (trate_abs >= rr1 .and. trate_abs < rr2) then
              cc1 = chron_ref(i,kc,2); cc2 = chron_ref(i,kc+1,2)
              t_closure = cc1 + (trate_abs-rr1)*(cc2-cc1)/(rr2-rr1)
            end if
          end do
        end if
      else
        t_closure = chron_ref(i,1,2)
      end if
      ! for the marker which already closure
      if (t_temp - t_closure .le. sterr .and. i_closure_if .eq. 1) then
        mark_chron_temp(i, kk) = t_closure
        cycle
      else if (t_closure - t_temp .le. sterr .and. i_closure_if .eq. 0)then
        mark_chron_temp(i, kk) = t_closure
        mark_chron_time(i, kk) = time
        cycle
      end if

      ! for the rock that still too hot to preserve the track
      if (t_temp - t_closure .gt. sterr) then
        i_closure_if = 0
        mark_chron_time(i, kk) = time
      ! for the rock that have been heated and cooled down,
      else if (t_closure - t_temp .gt. sterr) then
        i_closure_if = 1
        mark_chron_time(i, kk) = time
      end if

      mark_chron_if(i, kk) = i_closure_if
      mark_chron_temp(i, kk) = t_closure
    end do
!$OMP end parallel do

end do

return
end subroutine fl_thermchron

