! -*- F90 -*-
module phases
!! Phases the got special treatment in the code

integer, parameter :: koceandry = 10 ! dry oceanic crust
integer, parameter :: kcont1 = 2
integer, parameter :: kmant1 = 4
integer, parameter :: kocean1 = 3
integer, parameter :: kschist = 5
integer, parameter :: kcont2 = 6
integer, parameter :: kocean2 = 7
integer, parameter :: kmant2 = 8
integer, parameter :: kserp = 9
integer, parameter :: ksed1 = 1
integer, parameter :: ksed2 = 11
integer, parameter :: kweak = 12
integer, parameter :: keclg = 13
integer, parameter :: karc1 = 14
integer, parameter :: kweakmc = 15
integer, parameter :: khydmant = 16
integer, parameter :: kmetased = 17
integer, parameter :: kdrymant = 18
integer, parameter, dimension(5) :: mantle_phases = (/kmant1, kmant2, kserp, khydmant, kdrymant/)

! Orogeny lithologies (no change_phase rules attached -- intentional; see ARCHITECTURE §3.2/§3.3)
integer, parameter :: ksedrock = 20  ! sedimentary rock  (formerly borrowed 13 keclg)
integer, parameter :: kschisto = 21  ! schist            (formerly borrowed 14 karc1)
integer, parameter :: kultsch  = 22  ! ultra schist      (formerly borrowed 16 khydmant)
integer, parameter :: kdecoll  = 23  ! decollement       (formerly borrowed 18 kdrymant)
integer, parameter :: kbasemnt = 24  ! basement          (formerly borrowed 19)

end module phases
