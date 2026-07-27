subroutine fl_flexure
    use arrays
    use params
    implicit none

    double precision :: q(nx), Load(nx)
    double precision :: a_block(2, 2, nx), b_block(2, 2, nx), c_block(2, 2, nx)
    double precision :: r_block(2, nx), u_block(2, nx)
    double precision :: x(nx), dx(nx)
    double precision :: D(nx), E, nu, Te, rho_m, lambda, mu
    double precision :: cord_old_z, w_winkler_1, w_winkler_nx
    integer :: i, j, iph, bc_idx
    double precision :: a_i, c_i, b_i
    ! width (in nodes beyond the wall) of the band remesh re-stamps with
    ! pristine incoming stratigraphy; must match idist in remesh.f90
    integer, parameter :: nqref = 3

    ! 1. Scan for the flexural boundary condition (nbc = 200) on side 2 (bottom)
    bc_idx = 0
    do i = 1, nofbc
        if (nbc(i) .eq. 200 .and. nofside(i) .eq. 2) then
            bc_idx = i
            exit
        endif
    enddo

    if (bc_idx .eq. 0) return  ! Flexural boundary condition not active

    Te = bca(bc_idx)       ! Elastic thickness in meters
    rho_m = bcb(bc_idx)    ! Asthenosphere density in kg/m3

    ! Wrap the main computations in an OpenACC data region for stack arrays
    !$ACC data create(q, Load, a_block, b_block, c_block, r_block, u_block, x, dx, D)

    ! 2. Compute current vertical column weights (lithostatic loads)
    call get_column_weights(q)

    ! 3. If q_init is completely uninitialized (all zeros), initialize it with the starting column weights
    ! Handled cleanly on host and page-migrated automatically by CUDA Managed Memory
    if (all(q_init .eq. 0.d0)) then
        q_init = q
    endif

    ! Fix F (material-following load reference).  q_init is the column weight
    ! the plate is in equilibrium with; it is a property of the MATERIAL, not
    ! of the mesh column.  Lithosphere entering through a conveyor boundary
    ! arrives in isostatic equilibrium, so it must carry its own reference
    ! weight in with it.  Refresh the reference over the band that remesh
    ! re-stamps with pristine incoming stratigraphy (idist=2 -> elements
    ! 1..3 / nx-3..nx-1, i.e. nodes 1..1+nqref / nx-nqref..nx); rem_q_init
    ! then transports it inward on every rezoning, exactly like the phases.
    !
    ! Without this the initial LATERAL density contrast between the injected
    ! column and the interior columns is advected across the mesh and shows up
    ! as a travelling load.  Measured on taiwan_a g76_r000: incoming column
    ! 869.0 MPa vs interior 825.0 MPa = 44.0 MPa step at x=230 km, i.e.
    ! 1.33 km of spurious Airy subsidence propagating west at the convergence
    ! rate (trough observed at x = 230 - v*t to within 5 km at t=2/4/6/8).
    if (incoming_left .eq. 1) then
        do i = 1, min(1 + nqref, nx)
            q_init(i) = q(i)
        enddo
    endif
    if (incoming_right .eq. 1) then
        do i = max(nx - nqref, 1), nx
            q_init(i) = q(i)
        enddo
    endif

    ! Compute excess load
    !$ACC parallel loop async(1)
    !$OMP parallel do
    do i = 1, nx
        Load(i) = q(i) - q_init(i)
    enddo

    ! 4. Set up grid coordinates and spacings in X
    !$ACC parallel loop async(1)
    !$OMP parallel do
    do i = 1, nx
        x(i) = cord(nz, i, 1)
    enddo
    
    !$ACC parallel loop async(1)
    !$OMP parallel do
    do i = 1, nx-1
        dx(i) = x(i+1) - x(i)
        if (dx(i) .le. 0.d0) dx(i) = 1.d0  ! Safeguard against zero/negative spacing
    enddo

    ! 5. Compute local flexural rigidity D(i) based on Lamé parameters of the bottom elements
    !$ACC parallel loop private(iph, lambda, mu, E, nu) async(1)
    !$OMP parallel do private(iph, lambda, mu, E, nu)
    do i = 1, nx
        if (i .lt. nx) then
            iph = iphase(nz-1, i)
        else
            iph = iphase(nz-1, nx-1)
        endif

        ! Ensure phase is within valid range
        if (iph .lt. 1 .or. iph .gt. nphase) iph = 1

        lambda = rl(iph)
        mu = rm(iph)

        if (lambda + mu .gt. 0.d0) then
            E = mu * (3.0d0 * lambda + 2.0d0 * mu) / (lambda + mu)
            nu = lambda / (2.0d0 * (lambda + mu))
        else
            ! Fallback to standard mantle values if Lamé parameters are not set
            E = 1.0d11
            nu = 0.25d0
        endif

        if (nu .ge. 1.d0 .or. nu .le. -1.d0) nu = 0.25d0

        D(i) = (E * (Te**3)) / (12.0d0 * (1.0d0 - (nu**2)))
        if (D(i) .le. 0.d0) D(i) = 1.d15  ! Safeguard against zero/negative rigidity
    enddo

    ! 6. Construct the Block Tridiagonal System of Equations
    ! Boundary condition at left end: M_1 = 0, w_1 = q(1) / (rho_m * g)
    !$ACC serial async(1)
    a_block(:, :, 1) = 0.d0
    c_block(:, :, 1) = 0.d0
    b_block(1, 1, 1) = 1.d0
    b_block(1, 2, 1) = 0.d0
    b_block(2, 1, 1) = 0.d0
    b_block(2, 2, 1) = 1.d0
    
    r_block(1, 1) = 0.d0
    if (rho_m * g .gt. 0.d0) then
        w_winkler_1 = Load(1) / (rho_m * g)
    else
        w_winkler_1 = 0.d0
    endif
    r_block(2, 1) = w_winkler_1
    !$ACC end serial

    ! Inner nodes using non-uniform 2nd-order central differences
    !$ACC parallel loop private(a_i, c_i, b_i) async(1)
    !$OMP parallel do private(a_i, c_i, b_i)
    do i = 2, nx-1
        a_i = 2.d0 / (dx(i-1) * (dx(i) + dx(i-1)))
        c_i = 2.d0 / (dx(i) * (dx(i) + dx(i-1)))
        b_i = -(a_i + c_i)

        ! Diagonal block elements
        a_block(1, 1, i) = a_i
        a_block(1, 2, i) = 0.d0
        a_block(2, 1, i) = 0.d0
        a_block(2, 2, i) = a_i

        c_block(1, 1, i) = c_i
        c_block(1, 2, i) = 0.d0
        c_block(2, 1, i) = 0.d0
        c_block(2, 2, i) = c_i

        b_block(1, 1, i) = b_i
        b_block(1, 2, i) = -rho_m * g
        b_block(2, 1, i) = 1.0d0 / D(i)
        b_block(2, 2, i) = b_i

        r_block(1, i) = -Load(i)
        r_block(2, i) = 0.d0
    enddo

    ! Boundary condition at right end: M_nx = 0, w_nx = q(nx) / (rho_m * g)
    !$ACC serial async(1)
    a_block(:, :, nx) = 0.d0
    c_block(:, :, nx) = 0.d0
    b_block(1, 1, nx) = 1.d0
    b_block(1, 2, nx) = 0.d0
    b_block(2, 1, nx) = 0.d0
    b_block(2, 2, nx) = 1.d0
    
    r_block(1, nx) = 0.d0
    if (rho_m * g .gt. 0.d0) then
        w_winkler_nx = Load(nx) / (rho_m * g)
    else
        w_winkler_nx = 0.d0
    endif
    r_block(2, nx) = w_winkler_nx
    !$ACC end serial

    ! 7. Solve the block tridiagonal system directly on GPU
    !$ACC serial async(1)
    call solve_block_tridiagonal(nx, a_block, b_block, c_block, r_block, u_block)
    !$ACC end serial

    ! 8. Apply new deflections, update bottom boundary coordinates & velocities
    !$ACC parallel loop private(cord_old_z) async(1)
    !$OMP parallel do private(cord_old_z)
    do i = 1, nx
        cord_old_z = cord(nz, i, 2)
        cord(nz, i, 2) = zoriginal(nz, i) - u_block(2, i)
        vel(nz, i, 2) = (cord(nz, i, 2) - cord_old_z) / dt
    enddo

    !$ACC wait(1)
    !$ACC end data

end subroutine fl_flexure


subroutine rem_q_init
    ! Fix F, part 2: transport the flexural load reference across a rezoning.
    !
    ! Between remeshings the mesh is Lagrangian, so q_init(i) -- being index-
    ! attached -- travels with its material for free.  A remesh breaks that:
    ! every other field is interpolated from the old mesh onto the new one
    ! (rem_interpolate), but q_init was never in that list, so it stayed
    ! pinned to the node index while the material slid out from under it.
    ! In a conveyor run that silently freezes the reference in the Eulerian
    ! frame, which is what turns an advected lateral density contrast into a
    ! travelling load.  Interpolate it here, along the bottom row (q_init is
    ! a per-column quantity and fl_flexure indexes its grid by cord(nz,i,1)).
    use arrays
    use params
    implicit none
    double precision :: qold(nx), xo(nx), xn(nx), t
    integer :: i, k

    if (all(q_init .eq. 0.d0)) return   ! flexure not active yet

    qold = q_init
    do i = 1, nx
        xo(i) = cordo(nz, i, 1)
        xn(i) = cord(nz, i, 1)
    enddo

    k = 1
    do i = 1, nx
        if (xn(i) .le. xo(1)) then
            q_init(i) = qold(1)
        elseif (xn(i) .ge. xo(nx)) then
            q_init(i) = qold(nx)
        else
            do while (k .lt. nx-1 .and. xo(k+1) .lt. xn(i))
                k = k + 1
            enddo
            t = (xn(i) - xo(k)) / max(xo(k+1) - xo(k), 1.d-30)
            q_init(i) = qold(k) + t * (qold(k+1) - qold(k))
        endif
    enddo
end subroutine rem_q_init


subroutine flexure_report
    ! Prove the knob is connected: one line per output step giving the load
    ! anomaly the plate actually sees and the deflection it produced, plus
    ! how far the reference has drifted from its t=0 state.  Cheap, host-side,
    ! and written both to stdout and to flexure_diag.dat.
    use arrays
    use params
    implicit none
    double precision :: q(nx), Load(nx), w(nx)
    integer :: i, bc_idx, imax
    logical, save :: first = .true.

    bc_idx = 0
    do i = 1, nofbc
        if (nbc(i) .eq. 200 .and. nofside(i) .eq. 2) then
            bc_idx = i
            exit
        endif
    enddo
    if (bc_idx .eq. 0) return
    if (all(q_init .eq. 0.d0)) return

    call get_column_weights(q)
    do i = 1, nx
        Load(i) = q(i) - q_init(i)
        w(i) = zoriginal(nz, i) - cord(nz, i, 2)
    enddo

    imax = maxloc(abs(Load), 1)

    if (first) then
        open (77, file='flexure_diag.dat', status='unknown')
        write (77, '(A)') '# time_Myr  x_Loadmax_km  Load_max_MPa  ' // &
                          'x_wmax_km  w_max_m  Load(1)_MPa  Load(nx)_MPa'
        first = .false.
    else
        open (77, file='flexure_diag.dat', status='old', position='append')
    endif
    write (77, '(7(1PE14.6))') time/sec_year/1.d6, &
        cord(nz, imax, 1)/1.d3, Load(imax)/1.d6, &
        cord(nz, maxloc(w, 1), 1)/1.d3, maxval(w), &
        Load(1)/1.d6, Load(nx)/1.d6
    close (77)

    write (6, '(A,F7.3,A,F9.2,A,F8.1,A,F9.1,A,F8.1)') &
        ' [flexure] t=', time/sec_year/1.d6, ' Myr  Load_max=', Load(imax)/1.d6, &
        ' MPa @x=', cord(nz, imax, 1)/1.d3, ' km   w_max=', maxval(w), &
        ' m @x=', cord(nz, maxloc(w, 1), 1)/1.d3
end subroutine flexure_report


subroutine get_column_weights(q)
    use arrays
    use params
    implicit none
    double precision, intent(out) :: q(nx)
    double precision :: rho, dz, rho_b
    integer :: i, j
    double precision, external :: Eff_dens
    !$ACC routine(Eff_dens) seq

    !$ACC parallel loop private(rho, dz, rho_b, j) async(1)
    !$OMP parallel do private(rho, dz, rho_b, j)
    do i = 1, nx
        q(i) = 0.d0
        do j = 1, nz-1
            if (i .eq. 1) then
                rho = Eff_dens(j, 1)
            elseif (i .eq. nx) then
                rho = Eff_dens(j, nx-1)
            else
                rho = 0.5d0 * (Eff_dens(j, i-1) + Eff_dens(j, i))
            endif
            dz = cord(j, i, 2) - cord(j+1, i, 2)
            q(i) = q(i) + rho * g * dz
        enddo
        ! Fix A: remove the flexural-deflection load feedback.
        ! A downward deflection lengthens the meshed column, inflating q and
        ! cancelling the Winkler restoring force rho_m*g*w -> runaway subsidence
        ! (worst where the bottom density is high). Subtract the weight of the
        ! segment between the undeflected reference (zoriginal) and the current
        ! deflected bottom, using the bottom-element density, so Load carries
        ! only the tectonic anomaly.
        if (i .eq. 1) then
            rho_b = Eff_dens(nz-1, 1)
        elseif (i .eq. nx) then
            rho_b = Eff_dens(nz-1, nx-1)
        else
            rho_b = 0.5d0 * (Eff_dens(nz-1, i-1) + Eff_dens(nz-1, i))
        endif
        q(i) = q(i) - rho_b * g * (zoriginal(nz, i) - cord(nz, i, 2))
    enddo
end subroutine get_column_weights


subroutine solve_block_tridiagonal(n, a, b, c, r, u)
    !$ACC routine seq
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: a(2, 2, n), b(2, 2, n), c(2, 2, n), r(2, n)
    double precision, intent(out) :: u(2, n)
    
    double precision :: cp(2, 2, n), rp(2, n)
    double precision :: denom(2, 2), inv_denom(2, 2), det
    double precision :: tmp_v(2)
    integer :: i
    
    ! For i = 1: B_1 is identity, so B_1^-1 is identity.
    cp(:, :, 1) = c(:, :, 1)
    rp(:, 1) = r(:, 1)
    
    ! Forward sweep:
    do i = 2, n-1
        ! denom = B_i - A_i * cp_{i-1}
        denom(1, 1) = b(1, 1, i) - a(1, 1, i) * cp(1, 1, i-1)
        denom(1, 2) = b(1, 2, i) - a(1, 1, i) * cp(1, 2, i-1)
        denom(2, 1) = b(2, 1, i) - a(2, 2, i) * cp(2, 1, i-1)
        denom(2, 2) = b(2, 2, i) - a(2, 2, i) * cp(2, 2, i-1)
        
        ! Invert denom
        det = denom(1, 1) * denom(2, 2) - denom(1, 2) * denom(2, 1)
        if (abs(det) .lt. 1.d-30) det = sign(1.d-30, det)
        
        inv_denom(1, 1) = denom(2, 2) / det
        inv_denom(1, 2) = -denom(1, 2) / det
        inv_denom(2, 1) = -denom(2, 1) / det
        inv_denom(2, 2) = denom(1, 1) / det
        
        ! cp_i = inv_denom * c_i
        cp(1, 1, i) = inv_denom(1, 1) * c(1, 1, i)
        cp(1, 2, i) = inv_denom(1, 2) * c(2, 2, i)
        cp(2, 1, i) = inv_denom(2, 1) * c(1, 1, i)
        cp(2, 2, i) = inv_denom(2, 2) * c(2, 2, i)
        
        ! rp_i = inv_denom * (r_i - A_i * rp_{i-1})
        tmp_v(1) = r(1, i) - a(1, 1, i) * rp(1, i-1)
        tmp_v(2) = r(2, i) - a(2, 2, i) * rp(2, i-1)
        
        rp(1, i) = inv_denom(1, 1) * tmp_v(1) + inv_denom(1, 2) * tmp_v(2)
        rp(2, i) = inv_denom(2, 1) * tmp_v(1) + inv_denom(2, 2) * tmp_v(2)
    enddo
    
    ! For i = n: A_n is 0, denom = B_n = identity.
    u(:, n) = r(:, n)
    
    ! Backward substitution:
    do i = n-1, 1, -1
        u(1, i) = rp(1, i) - (cp(1, 1, i) * u(1, i+1) + cp(1, 2, i) * u(2, i+1))
        u(2, i) = rp(2, i) - (cp(2, 1, i) * u(1, i+1) + cp(2, 2, i) * u(2, i+1))
    enddo
end subroutine solve_block_tridiagonal
