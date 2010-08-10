subroutine init_marker

USE marker_data

use arrays
include 'precision.inc'
include 'params.inc'
include 'arrays.inc'

common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
real(8), dimension(9) :: x_tr,y_tr
parameter( min_elmarkers = 0, max_elmarkers = 12 )
parameter( onesixth = 0.1666666666666666666666)
parameter( fivesixth = 0.8333333333333333333333)

nphase_counter = 0
max_markers = 9*(nz-1)*(nx-1)

nelemts  = (nx-1)*(nz-1)
! define euler coordinate of the markers
! Distribute evenly first then randomize the distribution
! to start 9 markers per elements
nmarkers = 0
kk = 0
! zones with 9 markers per elements
! calculate the id (element number) of the zones of high res

!call random_seed
!write(333,*) 'Call to random_seed(), result may be stochastic'

do j = 1 , nz-1
    do i = 1 , nx-1
        dx = cord(j,i+1,1)-cord(j,i,1)
        dy = cord(j+1,i,2)-cord(j,i,2)

        x_tr(1) = cord(j,i,1) + dx*onesixth
        y_tr(1) = cord(j,i,2) + dy*onesixth
        x_tr(2) = cord(j,i,1) + dx*0.5
        y_tr(2) = cord(j,i,2) + dy*onesixth
        x_tr(3) = cord(j,i,1) + dx*fivesixth
        y_tr(3) = cord(j,i,2) + dy*onesixth

        x_tr(4) = cord(j,i,1) + dx*onesixth
        y_tr(4) = cord(j,i,2) + dy*0.5
        x_tr(5) = cord(j,i,1) + dx*0.5
        y_tr(5) = cord(j,i,2) + dy*0.5
        x_tr(6) = cord(j,i,1) + dx*fivesixth
        y_tr(6) = cord(j,i,2) + dy*0.5

        x_tr(7) = cord(j,i,1) + dx*onesixth
        y_tr(7) = cord(j,i,2) + dy*fivesixth
        x_tr(8) = cord(j,i,1) + dx*0.5
        y_tr(8) = cord(j,i,2) + dy*fivesixth
        x_tr(9) = cord(j,i,1) + dx*fivesixth
        y_tr(9) = cord(j,i,2) + dy*fivesixth


! randomize the new coordinates inside the element
        l = 1
        do while (l .le. 9)
            call random_number(rx)
            call random_number(ry)
            rx = 0.5 - rx
            ry = 0.5 - ry
            ddx = dx*rx/3
            ddy = dy*ry/3
            xx = x_tr(l)+ddx
            yy = y_tr(l)+ddy

            ii = i
            jj = j
            call euler2bar(xx,yy,bar1,bar2,ntr,ii,jj,inc)
            if(ntr.eq.0) cycle
      
! define the markers for each elements in columns (9*nelemts=nmarkers)
            l = l + 1
            kk = kk +1 

            mark(kk)%x = xx
            mark(kk)%y = yy
            mark(kk)%dead = 1 
            mark(kk)%ID = kk 
            mark(kk)%a1 = bar1
            mark(kk)%a2 = bar2
            mark(kk)%ntriag = ntr
            mark(kk)%phase = iphase(j,i)
            mark(kk)%maps   = aps(j,i)
            mark(kk)%meII   = 0. 
            mark(kk)%mpres  = 0.
            mark(kk)%mtemp  = 0.

            nphase_counter(j,i,mark(kk)%phase) = nphase_counter(j,i,mark(kk)%phase) + 1
        enddo
    enddo
enddo
nmarkers = kk
write(333,*) '# of markers', nmarkers
return
end subroutine init_marker

