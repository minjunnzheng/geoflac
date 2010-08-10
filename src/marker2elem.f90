subroutine marker2elem 
  USE marker_data

  use arrays
  include 'precision.inc'
  include 'params.inc'
  include 'arrays.inc'
  common /markers/ nfreemarkers,ndeadmarkers,xmpt(mnz*mnx*2,2,3)
  parameter( min_elmarkers = 0, max_elmarkers = 12 )
  integer kph(1)
  !type(marker):: mark(nmarkers)

  ! Interpolate marker properties into elements
  ! Find the triangle in which each marker belongs


  do j = 1 , nz-1
      do i = 1 , nx-1
          kinc = sum(nphase_counter(j,i,:))

          !  if there are too few markers in the element, create a new one
          if(kinc.le.4) call SysMsg('marker2elem: too few markers in the element, create a new one')
          do while (kinc.le.4)
              x1 = min(cord(j  ,i  ,1), cord(j+1,i  ,1))
              y1 = min(cord(j  ,i  ,2), cord(j  ,i+1,2))
              x2 = max(cord(j+1,i+1,1), cord(j  ,i+1,1))
              y2 = max(cord(j+1,i+1,2), cord(j+1,i  ,2))

              call random_number(rx)
              call random_number(ry)

              xx = x1 + (0.5-rx)*(x2-x1)
              yy = y1 + (0.5-ry)*(y2-y1)

              ! Calculate barycentic coordinates
              ii = i
              jj = j
              call euler2bar(xx,yy,bar1,bar2,ntr,ii,jj,inc) 

              if(inc.eq.0) cycle

              kinc = kinc + 1
              nmarkers = nmarkers + 1
              mark(nmarkers)%x = xx
              mark(nmarkers)%y = yy
              mark(nmarkers)%a1 = bar1
              mark(nmarkers)%a2 = bar2 
              mark(nmarkers)%maps = aps(j,i) 
              mark(nmarkers)%ntriag = ntr 
              mark(nmarkers)%dead = 1
              ! Calculate strain
              mark(nmarkers)%meII = strainII(j,i)
              ! Calculate pressure and temperature
              tmpr = 0.25*(temp(j,i)+temp(j+1,i)+temp(j,i+1)+temp(j+1,i+1))
              mark(nmarkers)%mpres = stressI(j,i)
              mark(nmarkers)%mtemp = tmpr 
              mark(nmarkers)%ID = nmarkers 
              ! assign phase to the new marker
              mark(nmarkers)%phase = iphase(j,i)
              nphase_counter(j,i,mark(nmarkers)%phase) = nphase_counter(j,i,mark(nmarkers)%phase) + 1
          enddo

          phase_ratio(j,i,1:nphase) = nphase_counter(j,i,1:nphase) / float(kinc)

          ! the phase of this element is the most abundant marker phase
          kph = maxloc(nphase_counter(j,i,:))
          iphase(j,i) = kph(1)

          ! sometimes there are more than one phases that are equally abundant
          maxphase = maxval(nphase_counter(j,i,:))
          nmax = count(nphase_counter(j,i,:) == maxphase)
          if(nmax .gt. 1) then
              write(*,*) 'elem has equally abundant marker phases:', i,j,nmax,nphase_counter(j,i,:)
              write(*,*) 'choosing the 1st maxloc as the phase'
          endif

      enddo
  enddo

  return
end subroutine marker2elem
