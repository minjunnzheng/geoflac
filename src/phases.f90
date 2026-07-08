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

! Orogeny lithologies 20-27 (no change_phase rules attached -- intentional).
! Full 12-19 -> 20-27 remap keeping the metamorphic pairs together (user scheme).
integer, parameter :: ksediment = 20  ! sediment           (was 12 kweak)
integer, parameter :: ksedrock  = 21  ! sedimentary rock   (was 13 keclg)
integer, parameter :: kmeta1a   = 22  ! metamorphic rock 1 (was 14 karc1)
integer, parameter :: kmeta1b   = 23  ! metamorphic rock 1 (was 15 kweakmc)
integer, parameter :: kmeta2a   = 24  ! metamorphic rock 2 (was 16 khydmant)
integer, parameter :: kmeta2b   = 25  ! metamorphic rock 2 (was 17 kmetased)
integer, parameter :: kdecoll   = 26  ! decollement        (was 18 kdrymant)
integer, parameter :: kbasemnt  = 27  ! basement           (was 19)

end module phases
