    subroutine re_simulation

    implicit none

    ! Variables
    real(8), parameter :: pi = 3.1415926535897932_8
    integer :: i=0,j=0,k=0,ll=0,mm=0

    logical, parameter :: init=.true.
    real(8), parameter :: Re=1.0d0, nu=0.002d0, t_end=2.0d0
    real(8) :: t_start
    !integer, parameter :: nx_file=256
    integer, parameter :: nx0=256, ny0=nx0, nz0=nx0, nxp0=nx0+2, nyp0=ny0+2, nzp0=nz0+2
    real(8), parameter :: lx=2.0d0*pi, ly=2.0d0*pi, lz=2.0d0*pi, dx=lx/nx0, dy=ly/ny0, dz=lz/nz0, dx2=dx*dx, dy2=dy*dy, dz2=dz*dz
    real(8), parameter :: xu(nx0+1)=[(i*dx, i=0, nx0)],          yu(ny0+2)=[((i+0.5)*dy, i=-1, ny0)],   zu(nz0+2)=[((i+0.5)*dz, i=-1, nz0)]
    real(8), parameter :: xv(nx0+2)=[((i+0.5d0)*dx, i=-1, nx0)], yv(ny0+1)=[(i*dy, i=0, ny0)],          zv(nz0+2)=[((i+0.5d0)*dz, i=-1, nz0)]
    real(8), parameter :: xw(nx0+2)=[((i+0.5d0)*dx, i=-1, nx0)], yw(ny0+2)=[((i+0.5d0)*dy, i=-1, ny0)], zw(nz0+1)=[(i*dz, i=0, nz0)]
    real(8), parameter :: xp(nx0+2)=[((i+0.5d0)*dx, i=-1, nx0)], yp(ny0+2)=[((i+0.5d0)*dy, i=-1, ny0)], zp(nz0+2)=[((i+0.5d0)*dz, i=-1, nz0)]

    integer, parameter :: nx=32, ny=nx, nz=nx, nxp=nx+2, nyp=ny+2, nzp=nz+2
    integer, parameter :: x0=16, y0=x0, z0=x0
    integer, parameter :: idx_xu(nx+1)=[(x0+i, i=0, nx)],   idx_yu(ny+2)=[(y0+i, i=0, ny+1)], idx_zu(nz+2)=[(z0+i, i=0, nz+1)]
    integer, parameter :: idx_xv(nx+2)=[(x0+i, i=0, nx+1)], idx_yv(ny+1)=[(y0+i, i=0, ny)],   idx_zv(nz+2)=[(z0+i, i=0, nz+1)]
    integer, parameter :: idx_xw(nx+2)=[(x0+i, i=0, nx+1)], idx_yw(ny+2)=[(y0+i, i=0, ny+1)], idx_zw(nz+1)=[(z0+i, i=0, nz)]
    integer, parameter :: idx_xp(nx+2)=[(x0+i, i=0, nx+1)], idx_yp(ny+2)=[(y0+i, i=0, ny+1)], idx_zp(nz+2)=[(z0+i, i=0, nz+1)]

    real(8) :: u(nx+1,ny+2,nz+2)=0,              v(nx+2,ny+1,nz+2)=0,              w(nx+2,ny+2,nz+1)=0,       p(nxp,nyp,nzp)=0
    real(8) :: u_star(nx+1,ny+2,nz+2)=0,         v_star(nx+2,ny+1,nz+2)=0,         w_star(nx+2,ny+2,nz+1)=0,  dp(nxp,nyp,nzp)=0
    real(8) :: rhs_x(nx+1,ny+2,nz+2)=0,          rhs_y(nx+2,ny+1,nz+2)=0,          rhs_z(nx+2,ny+2,nz+1)=0
    real(8) :: rhs_x_previous(nx+1,ny+2,nz+2)=0, rhs_y_previous(nx+2,ny+1,nz+2)=0, rhs_z_previous(nx+2,ny+2,nz+1)=0
    real(8) :: dpdx(nx+1,ny+2,nz+2)=0,           dpdy(nx+2,ny+1,nz+2)=0,           dpdz(nx+2,ny+2,nz+1)=0
    real(8) :: f_term_x=0,                       f_term_y=0,                       f_term_z=0
    !real(8) :: rhs_x_previous0(nx+1,ny+2,nz+2)=0,rhs_y_previous0(nx+2,ny+1,nz+2)=0,rhs_z_previous0(nx+2,ny+2,nz+1)=0
    !real(8) :: f_term_x(nx+1,ny+2,nz+2)=0,       f_term_y(nx+2,ny+1,nz+2)=0,       f_term_z(nx+2,ny+2,nz+1)=0

    real(8) :: dp_lu(nxp,nyp,nzp)=0, RHS_poisson0(nx*ny*nz)=0, dp_vec(nx*ny*nz)=0
    real(8), dimension (:,:,:), allocatable :: conv_x, conv_y, conv_z
    real(8), dimension (:,:,:), allocatable :: diff_x, diff_y, diff_z
    real(8) :: RHS_poisson_internal(nx,ny,nz), div(nx,ny,nz)

    !!!!!!!!!!!!!!! CN2-ADI scheme !!!!!!!!!!!!!!!
    real(8), dimension (:), allocatable :: A_low, A_d, A_up, A_up2, B_low, B_d, B_up, B_up2, C_low, C_d, C_up, C_up2 !CN2 scheme
    real(8), dimension (:), allocatable :: D_low, D_d, D_up, D_up2, E_low, E_d, E_up, E_up2, F_low, F_d, F_up, F_up2
    !real(8), dimension (:,:), allocatable :: A_mat, B_mat, C_mat, D_mat, E_mat, F_mat
    integer, dimension (:), allocatable :: A_ipiv, B_ipiv, C_ipiv, D_ipiv, E_ipiv, F_ipiv

    !!!!!!!!!!!!!!! big-simulation data !!!!!!!!!!!!!!!
    real(8) :: u_sub(nx+1,ny+2,nz+2)=0,          v_sub(nx+2,ny+1,nz+2)=0,          w_sub(nx+2,ny+2,nz+1)=0,       p_sub(nx+2,ny+2,nz+2)=0
    real(8) :: u_star_sub(nx+1,ny+2,nz+2)=0,     v_star_sub(nx+2,ny+1,nz+2)=0,     w_star_sub(nx+2,ny+2,nz+1)=0,  dp_sub(nx+2,ny+2,nz+2)=0
    real(8) :: RHS_poisson_sub(nx+2,ny+2,nz+2)=0,  p_sub_pre(nx+2,ny+2,nz+2)=0

    !!!!!!!!!!!!!!! boundary conditions !!!!!!!!!!!!!!!
    real(8) :: bx_u_1(ny+2,nz+2)=0, bx_u_nx(ny+2,nz+2)=0, by_u_1(nx+1,nz+2)=0, by_u_ny(nx+1,nz+2)=0, bz_u_1(nx+1,ny+2)=0, bz_u_nz(nx+1,ny+2)=0
    real(8) :: bx_v_1(ny+1,nz+2)=0, bx_v_nx(ny+1,nz+2)=0, by_v_1(nx+2,nz+2)=0, by_v_ny(nx+2,nz+2)=0, bz_v_1(nx+2,ny+1)=0, bz_v_nz(nx+2,ny+1)=0
    real(8) :: bx_w_1(ny+2,nz+1)=0, bx_w_nx(ny+2,nz+1)=0, by_w_1(nx+2,nz+1)=0, by_w_ny(nx+2,nz+1)=0, bz_w_1(nx+2,ny+2)=0, bz_w_nz(nx+2,ny+2)=0
    real(8) :: bx_p_1(nyp,nzp)=0,   bx_p_nx(nyp,nzp)=0,   by_p_1(nxp,nzp)=0,   by_p_ny(nxp,nzp)=0,   bz_p_1(nxp,nyp)=0,   bz_p_nz(nxp,nyp)=0

    real(8) :: tGet
    integer :: time_length, t_step, t_step_offset, plot_step=20, slice=nz/2+1
    real(8), dimension (:), allocatable :: time_array
    integer, dimension (:), allocatable :: x_range0, y_range0, z_range0
    integer :: x_range1(nx)=[(i, i=2, nx+1)], y_range1(ny)=[(i, i=2, ny+1)], z_range1(nz)=[(i, i=2, nz+1)]
    type(csr), allocatable :: LHS_poisson
    !type(coo), allocatable :: LHS_poisson_coo
    !real(8), dimension (:,:), allocatable :: LHS_poisson_den

    !!!!!!!!!!!!!!! temperary variables !!!!!!!!!!!!!!!
    !integer(8) :: sizeof_record, sizeof_record_sub, tempi1, tempi2
    !real(8) :: temp01, temp02, temp03!, temp04, temp05, temp06
    !real(8), dimension (:), allocatable :: temp11!, temp12, temp13, temp14, temp15, temp16
    real(8), dimension (:,:), allocatable :: temp21!, temp22, temp23, temp24, temp25, temp26
    real(8), dimension (:,:,:), allocatable :: temp31!, temp32, temp33!, temp34, temp35, temp36

    !!!!!!!!!!!!!!! INTEL mkl_pardiso !!!!!!!!!!!!!!!
    type(MKL_PARDISO_HANDLE) pt(64)
    integer :: maxfct=1, mnum=1, mtype=11, phase=13, n=nx*ny*nz, idum(nx*ny*nz), nrhs=1, iparm(64)=0, msglvl=0, error=0

    !!!!!!!!!!!!!!! INTEL mkl_dft !!!!!!!!!!!!!!!
    integer :: cstrides(4)=0, rstrides(4)=0
    type(DFTI_DESCRIPTOR), POINTER :: hand_f, hand_b
    !real(8) :: dft_out_r(nx,ny,nz), poisson_eigv(nx,ny,nz)=0
    !complex(8) :: dft_in_c(nx,ny,nz), dft_out_c(nx,ny,nz)
    real(8) :: dft_out_r1(nx,ny,nz), poisson_eigv(INT(nx/2.0)+1,ny,nz)=0
    complex(8) :: dft_out_c1(INT(nx/2.0)+1,ny,nz)
    integer :: status

    !!!!!!!!!!!!!!! INTEL mkl_tt !!!!!!!!!!!!!!!
    integer :: ipar_x(128), ipar_y(128), ipar_z(128)
    real(8) :: dpar_x(5*nx/2+2), dpar_y(5*nx/2+2), dpar_z(5*nz/2+2), rhs_tt(nx+1, ny+1, nz+1), eig_tt(nx, ny, nz)
    type(dfti_descriptor), pointer :: handle_x, handle_y, handle_z

    !!!!!!!!!!!!!!! HDF5 !!!!!!!!!!!!!!!
    character(len=1000) :: string_var, string_var2, big_DNS_file
    integer(8) :: h5f_whole, h5f_sub, h5f_slice, h5f_bc_interp, h5g_sub, h5g_slice, h5g_bc_interp

    !!!!!!!!!!!!!!! system_clock !!!!!!!!!!!!!!!
    REAL(8) :: system_clock_rate
    INTEGER :: c01,c02,c1,c2,cr,cm

    !!!!!!!!!!!!!!! re-simulation parameters !!!!!!!!!!!!!!!
    real(8) :: dt0=4e-3  !4.0d-3, 2.0d-3, 1.0d-3, 5.0d-4, 2.5d-4
    character(*), parameter :: timescheme="AB2"
    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall); pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    integer, parameter :: dp_flag=1, bc_x=2, bc_y=bc_x, bc_z=bc_x, pbc_x=3, pbc_y=3, pbc_z=3, sub_tstep=1
    logical, parameter :: using_Ustar=.false., TOffset=.true., restart=.true., BC_interp=.false.
    real(8), parameter :: noise=0
    real(8), dimension (:,:), allocatable :: err_vel, err_grad, rms_vel, rms_grad
    logical, parameter :: save_output=.false., LU_poisson=(.false. .and. nxp*nyp*nzp<=34**3), FFT_poisson=(pbc_x==1 .and. pbc_y==1 .and. pbc_z==1), DCT_poisson=(pbc_x/=1 .and. pbc_y/=1 .and. pbc_z/=1)

    call OMP_set_dynamic(.true.)
    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)
    call h5open_f(status)

    if (.not. TOffset) then
        t_start=0.0d0
    else
        t_start=1.0d0
    end if
    t_step_offset=nint(t_start/dt0)
    call get_command_argument (1, string_var, status=status)
    if (status==0) then
        read(string_var,*) dt0
    end if
    time_length=nint((t_end-t_start)/(dt0))
    allocate( time_array(0:time_length), err_vel(8,0:time_length), err_grad(24,0:time_length), rms_vel(4,0:time_length), rms_grad(12,0:time_length) )
    time_array=[(0.0d0+dt0*i, i=t_step_offset, time_length+t_step_offset)]
    allocate(temp31(nx+2,ny+2,nz+2))

    write (string_var2,'(A, "_", ES5.0E1)') trim(timescheme), dt0
    if (.not. dp_flag) string_var2=trim(string_var2)//"_p"
    if (pbc_x==2) string_var2=trim(string_var2)//"_D"
    if (pbc_x==3) string_var2=trim(string_var2)//"_N"
    if (pbc_x==4) string_var2=trim(string_var2)//"_Dg"
    if (using_Ustar) string_var2=trim(string_var2)//"_Ustar"
    if (.not. using_Ustar) string_var2=trim(string_var2)//"_U"
    if (TOffset) string_var2=trim(string_var2)//"_TOffset"
    if (restart) string_var2=trim(string_var2)//"_restart"
    if (noise/=0.0d0) write (string_var2,'(A, "_noise", ES6.0E2)') trim(string_var2), noise
    if (sub_tstep/=1) write (string_var2,'(A, "_sub_tstep", I0)') trim(string_var2), sub_tstep
    if (BC_interp) string_var2=trim(string_var2)//"_BCinterp"

    if (bc_x==1) then
        !allocate( x_range0(nx) )
        x_range0=[ (i, i=1, nx) ]
    else
        !allocate( x_range0(nx-1) )
        x_range0=[ (i, i=2, nx) ]
    end if
    if (bc_y==1) then
        !allocate( y_range0(ny) )
        y_range0=[ (i, i=1, ny) ]
    else
        !allocate( y_range0(ny-1) )
        y_range0=[ (i, i=2, ny) ]
    end if
    if (bc_z==1) then
        !allocate( z_range0(nz) )
        z_range0=[ (i, i=1, nz) ]
    else
        !allocate( z_range0(nx-1) )
        z_range0=[ (i, i=2, nz) ]
    end if

    if (bc_x==1) then
        allocate( conv_x(nx,ny,nz),   diff_x(nx,ny,nz)   )
    else
        allocate( conv_x(nx-1,ny,nz), diff_x(nx-1,ny,nz) )
    end if

    if (bc_y==1) then
        allocate( conv_y(nx,ny,nz),   diff_y(nx,ny,nz)   )
    else
        allocate( conv_y(nx,ny-1,nz), diff_y(nx,ny-1,nz) )
    end if

    if (bc_z==1) then
        allocate( conv_z(nx,ny,nz),   diff_z(nx,ny,nz)   )
    else
        allocate( conv_z(nx,ny,nz-1), diff_z(nx,ny,nz-1) )
    end if

    if (timescheme=="AB2-CN") then !assuming bc_x=bc_y=bc_z, nx=ny=nz
        if (allocated(A_low)) deallocate(A_low, A_d, A_up, A_up2, A_ipiv, B_low, B_d, B_up, B_up2, B_ipiv)
        if (allocated(C_low)) deallocate(C_low, C_d, C_up, C_up2, C_ipiv, D_low, D_d, D_up, D_up2, D_ipiv)
        if (allocated(E_low)) deallocate(E_low, E_d, E_up, E_up2, E_ipiv, F_low, F_d, F_up, F_up2, F_ipiv)
        allocate(A_low(nx), A_d(nx+1), A_up(nx), A_up2(nx-1), A_ipiv(nx+1), B_low(nx+1), B_d(nx+2), B_up(nx+1), B_up2(nx), B_ipiv(nx+2))
        allocate(C_low(ny), C_d(ny+1), C_up(ny), C_up2(ny-1), C_ipiv(ny+1), D_low(ny+1), D_d(ny+2), D_up(ny+1), D_up2(ny), D_ipiv(ny+2))
        allocate(E_low(nz), E_d(nz+1), E_up(nz), E_up2(nz-1), E_ipiv(nz+1), F_low(nz+1), F_d(nz+2), F_up(nz+1), F_up2(nz), F_ipiv(nz+2))

        call mat_CN_tridiagonal(nx, bc_x, dx2, dx, A_low, A_d, A_up, B_low, B_d, B_up, dt0, nu)
        call mat_CN_tridiagonal(ny, bc_y, dy2, dy, C_low, C_d, C_up, D_low, D_d, D_up, dt0, nu)
        call mat_CN_tridiagonal(nz, bc_z, dz2, dz, E_low, E_d, E_up, F_low, F_d, F_up, dt0, nu)

        !CN matrix should be a diagonally dominant tridiagonal coefficient matrix (?)
        call dttrfb( A_low, A_d, A_up )
        call dttrfb( B_low, B_d, B_up )
        call dttrfb( C_low, C_d, C_up )
        call dttrfb( D_low, D_d, D_up )
        call dttrfb( E_low, E_d, E_up )
        call dttrfb( F_low, F_d, F_up )
    end if

    !!!! LU-based FD poisson solver
    if (LU_poisson) then
        LHS_poisson=Poisson_LHS_staggered(nx, ny, nz, dx2, dy2, dz2, pbc_x, pbc_y, pbc_z, dx, dy, dz)
        !LHS_poisson_coo=LHS_poisson_coo%from_csr(LHS_poisson)
        !LHS_poisson_den=LHS_poisson_coo%to_den()

        !OPEN(10, file="poisson_eq.dat", form="unformatted")
        !WRITE(10) [LHS_poisson_den]
        !CLOSE(10)

        DO i = 1, 64
            pt(i)%DUMMY = 0
        END DO
        phase=12
        iparm(1) = 0 ! no solver default
        call pardisoinit (pt, mtype, iparm)
        print *,iparm(1)
        iparm(2) = 3 ! parallel (OpenMP) version of the nested dissection algorithm
        !iparm(5) = 2 ! no user fill-in reducing permutation
        !iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        !iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
        !iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(25)=1 ! Intel MKL PARDISO uses a parallel algorithm for the solve step.

        print *, "**************************************"
        print *, "LU of LHS_poisson start..."
        CALL SYSTEM_CLOCK(c1)
        call pardiso (pt, maxfct, mnum, mtype, phase, n, LHS_poisson%value, LHS_poisson%ia, LHS_poisson%ja, &
            idum, nrhs, iparm, msglvl, RHS_poisson0, dp_vec, error)
        CALL SYSTEM_CLOCK(c2)
        print '(" LU of LHS_poisson completed: ", F8.4, " second")', (c2-c1)/system_clock_rate
        print *, "**************************************"
    end if

    !!!! FFT-based FD poisson solver
    if (FFT_poisson) then
        ! Configure Forward Descriptor
        hand_f => null()
        hand_b => null()
        cstrides = [0, 1, INT(nx/2.0)+1, ny*(INT(nx/2.0)+1)]
        rstrides = [0, 1, nx,     ny*nx]
        !3D MKL dft = fft_z( fft_y (fft_x) ), same in MATLAB
        print *,"Configure DFTI descriptor for fordward transform"
        !status = DftiCreateDescriptor(hand_f, DFTI_DOUBLE, DFTI_COMPLEX, 3, [nx,ny,nz])
        status = DftiCreateDescriptor(hand_f, DFTI_DOUBLE, DFTI_REAL, 3, [nx,ny,nz])
        status = DftiSetValue(hand_f, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        status = DftiSetValue(hand_f, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
        status = DftiSetValue(hand_f, DFTI_INPUT_STRIDES, rstrides)
        status = DftiSetValue(hand_f, DFTI_OUTPUT_STRIDES, cstrides)
        status = DftiSetValue(hand_f, DFTI_BACKWARD_SCALE, 1.0d0/(nx*ny*nz))
        status = DftiCommitDescriptor(hand_f)

        print *,"Configure DFTI descriptor for backward transform"
        !status = DftiCreateDescriptor(hand_b, DFTI_DOUBLE, DFTI_COMPLEX, 3, [nx,ny,nz])
        status = DftiCreateDescriptor(hand_b, DFTI_DOUBLE, DFTI_REAL, 3, [nx,ny,nz])
        status = DftiSetValue(hand_b, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        status = DftiSetValue(hand_b, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
        status = DftiSetValue(hand_b, DFTI_INPUT_STRIDES, cstrides)
        status = DftiSetValue(hand_b, DFTI_OUTPUT_STRIDES, rstrides)
        status = DftiSetValue(hand_b, DFTI_BACKWARD_SCALE, 1.0d0/(nx*ny*nz))
        status = DftiCommitDescriptor(hand_b)

        do i=1,nx/2+1
            poisson_eigv(i,:,:)=(sin(pi*(i-1)/nx)/dx)**2
        end do

        do i=1,ny
            poisson_eigv(:,i,:)=poisson_eigv(:,i,:)+(sin(pi*(i-1)/ny)/dy)**2
        end do

        do i=1,nz
            poisson_eigv(:,:,i)=poisson_eigv(:,:,i)+(sin(pi*(i-1)/nz)/dz)**2
        end do
        !poisson_eigv(1,1,1)=1.0d0
        poisson_eigv=-4.0d0*poisson_eigv
        poisson_eigv(1,1,1) = IEEE_VALUE (1.0d0, IEEE_POSITIVE_INF)
    end if

    !!!! DCT-based FD poisson solver (pbc=1 and 4 are not implemented and tested)
    if (DCT_poisson) then
        if (pbc_x==4) then
            call d_init_trig_transform(nx, MKL_SINE_TRANSFORM, ipar_x, dpar_x, status)
        else
            call d_init_trig_transform(nx, MKL_STAGGERED_COSINE_TRANSFORM, ipar_x, dpar_x, status)
        end if
        if (pbc_y==4) then
            call d_init_trig_transform(ny, MKL_SINE_TRANSFORM, ipar_y, dpar_y, status)
        else
            call d_init_trig_transform(ny, MKL_STAGGERED_COSINE_TRANSFORM, ipar_y, dpar_y, status)
        end if
        if (pbc_z==4) then
            call d_init_trig_transform(nz, MKL_SINE_TRANSFORM, ipar_z, dpar_z, status)
        else
            call d_init_trig_transform(nz, MKL_STAGGERED_COSINE_TRANSFORM, ipar_z, dpar_z, status)
        end if
        call d_commit_trig_transform(rhs_tt(:,1,1),handle_x,ipar_x,dpar_x,status)
        call d_commit_trig_transform(rhs_tt(1,:,1),handle_y,ipar_y,dpar_y,status)
        call d_commit_trig_transform(rhs_tt(1,1,:),handle_z,ipar_z,dpar_z,status)

        do i=1,nx
            if (pbc_x==2) then
                eig_tt(i,:,:)=(sin( i   *pi/2/nx))**2/dx2
            elseif (pbc_x==3) then
                eig_tt(i,:,:)=(sin((i-1)*pi/2/nx))**2/dx2
            elseif (pbc_x==4) then
                eig_tt(i,:,:)=(sin(i*pi/2/(nx+1)))**2/dx2
            end if
        end do

        do i=1,ny
            if (pbc_y==2) then
                eig_tt(:,i,:)=eig_tt(:,i,:)+(sin( i   *pi/2/ny))**2/dy2
            elseif (pbc_y==3) then
                eig_tt(:,i,:)=eig_tt(:,i,:)+(sin((i-1)*pi/2/ny))**2/dy2
            elseif (pbc_y==4) then
                eig_tt(:,i,:)=eig_tt(:,i,:)+(sin(i*pi/2/(ny+1)))**2/dy2
            end if
        end do

        do i=1,nz
            if (pbc_z==2) then
                !eig_tt(:,:,i)=eig_tt(:,:,i)+(sin( i   *pi/2/nz))**2/dz2
                eig_tt(:,:,i)=eig_tt(:,:,i)+(sin((nz-i+1)*pi/2/nz))**2/dz2 !!!hard-coded to speed up used in eigenvalue normalization
            elseif (pbc_z==3) then
                eig_tt(:,:,i)=eig_tt(:,:,i)+(sin((i-1)*pi/2/nz))**2/dz2
            elseif (pbc_y==4) then
                eig_tt(:,:,i)=eig_tt(:,:,i)+(sin(i*pi/2/(nz+1)))**2/dz2
            end if
        end do

        !if (eig_tt(1,1,1)==0.0d0) eig_tt(1,1,1)=1.0d0
        eig_tt=-4.0d0*eig_tt
        if (eig_tt(1,1,1)==0.0d0) eig_tt(1,1,1) = IEEE_VALUE (1.0d0, IEEE_POSITIVE_INF)
    end if

    if (timescheme=="Euler") then
        if (dp_flag) then
            write (big_DNS_file,'(A, "_result.MARCC/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_dp_x0_", I0, "_nx0_", I0, "_sub.h5")') trim("AB2"), nx0, dt0, trim("AB2"), x0, nx
        else
            write (big_DNS_file,'(A, "_result.MARCC/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_p_x0_", I0, "_nx0_", I0, "_sub.h5")') trim("AB2"), nx0, dt0, trim("AB2"), x0, nx
        end if
    else
        if (dp_flag) then
            write (big_DNS_file,'(A, "_result.MARCC/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_dp_x0_", I0, "_nx0_", I0, "_sub.h5")') trim(timescheme), nx0, dt0, trim(timescheme), x0, nx
        else
            write (big_DNS_file,'(A, "_result.MARCC/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_p_x0_", I0, "_nx0_", I0, "_sub.h5")') trim(timescheme), nx0, dt0, trim(timescheme), x0, nx
        end if
    end if

    call h5fopen_f(big_DNS_file, H5F_ACC_RDONLY_F, h5f_sub, status)

    if (BC_interp) then
        if (using_Ustar==.true. .and. TOffset==.true. .and. restart==.true. .and. noise==0.0d0 .and. sub_tstep==1 .and. timescheme=="AB2" .and. dp_flag==1) then
            write (string_var,'(A, "_result.MARCC/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_dp_x0_", I0, "_nx0_", I0, "_sub_BC_spline.h5")') trim("AB2"), nx0, dt0, trim("AB2"), x0, nx
            call h5fopen_f(string_var, H5F_ACC_RDONLY_F, h5f_bc_interp, status)
        else
            print *, 'BC_interp errors!!! re_simulation.f90: Line 361!!!'
            stop(1)
        end if
    end if

    do t_step=0,time_length
        tGet=time_array(t_step)
        print *,''
        print '(" t_step ", I6, "        tGet ", F7.4)', t_step, tGet
        CALL SYSTEM_CLOCK(c01)

        if (t_step==0) then
            write (string_var,'("t_", F0.4)') tGet
            call h5gopen_f(h5f_sub, string_var, h5g_sub, status)
            call h5ltread_dataset_double_f(h5g_sub, 'u_sub', u_sub, [nx+1_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'v_sub', v_sub, [nx+2_8,ny+1_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'w_sub', w_sub, [nx+2_8,ny+2_8,nz+1_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'p_sub', p_sub, [nx+2_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'rhs_x_previous', rhs_x_previous, [nx+1_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'rhs_y_previous', rhs_y_previous, [nx+2_8,ny+1_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'rhs_z_previous', rhs_z_previous, [nx+2_8,ny+2_8,nz+1_8], status)
            !call h5ltread_dataset_double_f(h5g_sub, 'u_star_sub', u_star_sub, [nx+1_8,ny+2_8,nz+2_8], status)
            !call h5ltread_dataset_double_f(h5g_sub, 'v_star_sub', v_star_sub, [nx+2_8,ny+1_8,nz+2_8], status)
            !call h5ltread_dataset_double_f(h5g_sub, 'w_star_sub', w_star_sub, [nx+2_8,ny+2_8,nz+1_8], status)
            call h5gclose_f( h5g_sub, status)

            u=u_sub; v=v_sub; w=w_sub; p=p_sub;

            !call TGV(xu(idx_xu), yu(idx_yu), zu(idx_zu), 0.0d0, nu, u)
            !call TGV(xv(idx_xv), yv(idx_yv), zv(idx_zv), 0.0d0, nu, v=v)
            !call TGV(xp(idx_xp), yp(idx_yp), zp(idx_zp), 0.0d0, nu, p=p)
        else
            if (.not. dp_flag) then
                p_sub_pre=p_sub
            end if
            
            write (string_var,'("t_", F0.4)') tGet
            call h5gopen_f(h5f_sub, string_var, h5g_sub, status)

            call h5ltread_dataset_double_f(h5g_sub, 'u_sub', u_sub, [nx+1_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'v_sub', v_sub, [nx+2_8,ny+1_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'w_sub', w_sub, [nx+2_8,ny+2_8,nz+1_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'p_sub', p_sub, [nx+2_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'u_star_sub', u_star_sub, [nx+1_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'v_star_sub', v_star_sub, [nx+2_8,ny+1_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'w_star_sub', w_star_sub, [nx+2_8,ny+2_8,nz+1_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'dp_sub', dp_sub, [nx+2_8,ny+2_8,nz+2_8], status)
            call h5ltread_dataset_double_f(h5g_sub, 'RHS_poisson_sub', RHS_poisson_sub, [nx+2_8,ny+2_8,nz+2_8], status)
            if (timescheme=="AB2-CN") then
                call h5ltread_dataset_double_f(h5g_sub, 'bx_u_1s', bx_u_1, [ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bx_u_nxs', bx_u_nx, [ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bx_v_1s', bx_v_1, [ny+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bx_v_nxs', bx_v_nx, [ny+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bx_w_1s', bx_w_1, [ny+2_8,nz+1_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bx_w_nxs', bx_w_nx, [ny+2_8,nz+1_8], status)

                call h5ltread_dataset_double_f(h5g_sub, 'by_u_1s', by_u_1, [nx+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'by_u_nys', by_u_ny, [nx+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'by_v_1s', by_v_1, [nx+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'by_v_nys', by_v_ny, [nx+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'by_w_1s', by_w_1, [nx+2_8,nz+1_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'by_w_nys', by_w_ny, [nx+2_8,nz+1_8], status)

                call h5ltread_dataset_double_f(h5g_sub, 'bz_u_1s', bz_u_1, [nx+1_8,ny+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bz_u_nzs', bz_u_nz, [nx+1_8,ny+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bz_v_1s', bz_v_1, [nx+2_8,ny+1_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bz_v_nzs', bz_v_nz, [nx+2_8,ny+1_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bz_w_1s', bz_w_1, [nx+2_8,ny+2_8], status)
                call h5ltread_dataset_double_f(h5g_sub, 'bz_w_nzs', bz_w_nz, [nx+2_8,ny+2_8], status)
            end if

            call h5gclose_f( h5g_sub, status)

            if (BC_interp) then
                call h5gopen_f(h5f_bc_interp, string_var, h5g_bc_interp, status)
                call h5ltread_dataset_double_f(h5g_bc_interp, 'u_star_sub', u_star_sub, [nx+1_8,ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_bc_interp, 'v_star_sub', v_star_sub, [nx+2_8,ny+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5g_bc_interp, 'w_star_sub', w_star_sub, [nx+2_8,ny+2_8,nz+1_8], status)
                call h5ltread_dataset_double_f(h5g_bc_interp, 'dp_sub', dp_sub, [nx+2_8,ny+2_8,nz+2_8], status)
                call h5gclose_f( h5g_bc_interp, status)
            end if

            if (.not. dp_flag) then
                dp_sub=p_sub
            end if

            if (timescheme=="AB2-CN") then
                call get_pr_bc(dp_sub, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)
            else if (using_Ustar) then
                call get_velpr_bc(u_star_sub, v_star_sub, w_star_sub, dp_sub, bc_x, bc_y, bc_z, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, &
                    bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz, bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz, &
                    bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz, dt0, dp_flag, using_Ustar,p_sub_pre)
            else
                call get_velpr_bc(u_sub, v_sub, w_sub, dp_sub, bc_x, bc_y, bc_z, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, &
                    bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz, bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz, &
                    bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz, dt0, dp_flag, using_Ustar,p_sub_pre)
            end if

            if (noise>0 .and. t_step==1) then
                bx_u_1 =whiteNoise2(bx_u_1, noise)
                bx_u_nx=whiteNoise2(bx_u_nx,noise)
                bx_v_1 =whiteNoise2(bx_v_1, noise)
                bx_v_nx=whiteNoise2(bx_v_nx,noise)
                bx_w_1 =whiteNoise2(bx_w_1, noise)
                bx_w_nx=whiteNoise2(bx_w_nx,noise)

                by_u_1 =whiteNoise2(by_u_1, noise)
                by_u_ny=whiteNoise2(by_u_ny,noise)
                by_v_1 =whiteNoise2(by_v_1, noise)
                by_v_ny=whiteNoise2(by_v_ny,noise)
                by_w_1 =whiteNoise2(by_w_1, noise)
                by_w_ny=whiteNoise2(by_w_ny,noise)

                bz_u_1 =whiteNoise2(bz_u_1, noise)
                bz_u_nz=whiteNoise2(bz_u_nz,noise)
                bz_v_1 =whiteNoise2(bz_v_1, noise)
                bz_v_nz=whiteNoise2(bz_v_nz,noise)
                bz_w_1 =whiteNoise2(bz_w_1, noise)
                bz_w_nz=whiteNoise2(bz_w_nz,noise)

                !bx_p_1 =whiteNoise2(bx_p_1, noise)
                !bx_p_nx=whiteNoise2(bx_p_nx,noise)
                !by_p_1 =whiteNoise2(by_p_1, noise)
                !by_p_ny=whiteNoise2(by_p_ny,noise)
                !bz_p_1 =whiteNoise2(bz_p_1, noise)
                !bz_p_nz=whiteNoise2(bz_p_nz,noise)
            end if

            !!! Convection
            !$omp parallel sections if (nx<40)
            !$omp section
            call cal_conv(u, v, w, bc_x, bc_y, bc_z, dx, dy, dz, conv_x, conv_y, conv_z)

            !!! Diffusion
            !$omp section
            call cal_diff(u, v, w, bc_x, bc_y, bc_z, dx2, dy2, dz2, diff_x, diff_y, diff_z)

            !!! Pressure gradients
            !$omp section
            if (dp_flag) then
                dpdx=diff(p,1,1)/dx;
                dpdy=diff(p,1,2)/dy;
                dpdz=diff(p,1,3)/dz;
            end if

            f_term_x=0; f_term_y=0; f_term_z=0;
            !$omp end parallel sections

            !!! Time-advancement
            if (timescheme=="Euler") then
                ! diff terms + conv terms
                !$omp parallel sections if (nx<40)
                !$omp section
                rhs_x(x_range0,y_range1,z_range1)=nu*diff_x-conv_x;
                u_star=dt0*(1.0d0*rhs_x-dp_flag*dpdx+1.0d0*f_term_x)+u;
                !rhs_x_previous=0;
                !$omp section
                rhs_y(x_range1,y_range0,z_range1)=nu*diff_y-conv_y;
                v_star=dt0*(1.0d0*rhs_y-dp_flag*dpdy+1.0d0*f_term_y)+v;
                !rhs_y_previous=0;
                !$omp section
                rhs_z(x_range1,y_range1,z_range0)=nu*diff_z-conv_z;
                w_star=dt0*(1.0d0*rhs_z-dp_flag*dpdz+1.0d0*f_term_z)+w;
                !rhs_z_previous=0;
                !$omp end parallel sections

                call vel_bc_staggered(u_star,v_star,w_star,&
                    bx_u_1,bx_u_nx,by_u_1,by_u_ny,bz_u_1,bz_u_nz,&
                    bx_v_1,bx_v_nx,by_v_1,by_v_ny,bz_v_1,bz_v_nz,&
                    bx_w_1,bx_w_nx,by_w_1,by_w_ny,bz_w_1,bz_w_nz,&
                    bc_x,bc_y,bc_z);
            else if (timescheme=="AB2") then
                ! diff terms + conv terms

                ! prediction
                if (.not. restart .and. t_step==1) then
                    !$omp parallel sections if (nx<40)
                    !$omp section
                    rhs_x(x_range0,y_range1,z_range1)=nu*diff_x-conv_x;
                    u_star=dt0*(1.0d0*rhs_x-dp_flag*dpdx+1.0d0*f_term_x)+u;
                    rhs_x_previous=rhs_x;
                    !$omp section
                    rhs_y(x_range1,y_range0,z_range1)=nu*diff_y-conv_y;
                    v_star=dt0*(1.0d0*rhs_y-dp_flag*dpdy+1.0d0*f_term_y)+v;
                    rhs_y_previous=rhs_y;
                    !$omp section
                    rhs_z(x_range1,y_range1,z_range0)=nu*diff_z-conv_z;
                    w_star=dt0*(1.0d0*rhs_z-dp_flag*dpdz+1.0d0*f_term_z)+w;
                    rhs_z_previous=rhs_z;
                    !$omp end parallel sections
                else
                    !$omp parallel sections if (nx<40)
                    !$omp section
                    rhs_x(x_range0,y_range1,z_range1)=nu*diff_x-conv_x;
                    u_star=dt0*(1.5d0*rhs_x-0.5d0*rhs_x_previous-dp_flag*dpdx+1.0d0*f_term_x)+u;
                    rhs_x_previous=rhs_x;
                    !$omp section
                    rhs_y(x_range1,y_range0,z_range1)=nu*diff_y-conv_y;
                    v_star=dt0*(1.5d0*rhs_y-0.5d0*rhs_y_previous-dp_flag*dpdy+1.0d0*f_term_y)+v;
                    rhs_y_previous=rhs_y;
                    !$omp section
                    rhs_z(x_range1,y_range1,z_range0)=nu*diff_z-conv_z;
                    w_star=dt0*(1.5d0*rhs_z-0.5d0*rhs_z_previous-dp_flag*dpdz+1.0d0*f_term_z)+w;
                    rhs_z_previous=rhs_z;
                    !$omp end parallel sections

                end if

                call vel_bc_staggered(u_star,v_star,w_star,&
                    bx_u_1,bx_u_nx,by_u_1,by_u_ny,bz_u_1,bz_u_nz,&
                    bx_v_1,bx_v_nx,by_v_1,by_v_ny,bz_v_1,bz_v_nz,&
                    bx_w_1,bx_w_nx,by_w_1,by_w_ny,bz_w_1,bz_w_nz,&
                    bc_x,bc_y,bc_z);
            else if (timescheme=="AB2-CN") then
                ! prediction
                if ( .not. restart .and. t_step==1 ) then
                    ! diff terms + conv terms
                    !$omp parallel sections if (nx<40)
                    !$omp section
                    rhs_x(x_range0,y_range1,z_range1)=0.5d0*nu*diff_x-conv_x;
                    u_star=dt0*(rhs_x-dp_flag*dpdx+1.0d0*f_term_x)+u;
                    rhs_x_previous(x_range0,y_range1,z_range1)=-conv_x;
                    !$omp section
                    rhs_y(x_range1,y_range0,z_range1)=0.5d0*nu*diff_y-conv_y;
                    v_star=dt0*(rhs_y-dp_flag*dpdy+1.0d0*f_term_y)+v;
                    rhs_y_previous(x_range1,y_range0,z_range1)=-conv_y;
                    !$omp section
                    rhs_z(x_range1,y_range1,z_range0)=0.5d0*nu*diff_z-conv_z;
                    w_star=dt0*(rhs_z-dp_flag*dpdz+1.0d0*f_term_z)+w;
                    rhs_z_previous(x_range1,y_range1,z_range0)=-conv_z;
                    !$omp end parallel sections
                else
                    ! diff terms + conv terms
                    !$omp parallel sections if (nx<40)
                    !$omp section
                    rhs_x(x_range0,y_range1,z_range1)=0.5d0*nu*diff_x-1.5d0*conv_x;
                    u_star=dt0*(rhs_x-0.5d0*rhs_x_previous-dp_flag*dpdx+1.0d0*f_term_x)+u;
                    rhs_x_previous(x_range0,y_range1,z_range1)=-conv_x;
                    !$omp section
                    rhs_y(x_range1,y_range0,z_range1)=0.5d0*nu*diff_y-1.5d0*conv_y;
                    v_star=dt0*(rhs_y-0.5d0*rhs_y_previous-dp_flag*dpdy+1.0d0*f_term_y)+v;
                    rhs_y_previous(x_range1,y_range0,z_range1)=-conv_y;
                    !$omp section
                    rhs_z(x_range1,y_range1,z_range0)=0.5d0*nu*diff_z-1.5d0*conv_z;
                    w_star=dt0*(rhs_z-0.5d0*rhs_z_previous-dp_flag*dpdz+1.0d0*f_term_z)+w;
                    rhs_z_previous(x_range1,y_range1,z_range0)=-conv_z;
                    !$omp end parallel sections
                end if

                call vel_bc_staggered_CN2(u_star, v_star, w_star, bx_u_1, bx_u_nx, bx_v_1, bx_v_nx, bx_w_1, bx_w_nx, bc_x, 1)
                !$omp parallel do if (nx<40)
                do k=1,nz+2
                    call dttrsb( A_low, A_d, A_up, u_star(:,:,k) )
                    call dttrsb( B_low, B_d, B_up, v_star(:,:,k) )
                    if (k<=nz+1) call dttrsb( B_low, B_d, B_up, w_star(:,:,k) )
                end do
                !$omp end parallel do

                call vel_bc_staggered_CN2(u_star, v_star, w_star, by_u_1, by_u_ny, by_v_1, by_v_ny, by_w_1, by_w_ny, bc_y, 2)
                !$omp parallel do private(temp21) if (nx<40)
                do k=1,nz+2
                    temp21=transpose(u_star(:,:,k))
                    call dttrsb( D_low, D_d, D_up, temp21 )
                    u_star(:,:,k)=transpose(temp21)

                    temp21=transpose(v_star(:,:,k))
                    call dttrsb( C_low, C_d, C_up, temp21 )
                    v_star(:,:,k)=transpose(temp21)

                    if (k<=nz+1) then
                        temp21=transpose(w_star(:,:,k))
                        call dttrsb( D_low, D_d, D_up, temp21 )
                        w_star(:,:,k)=transpose(temp21)
                    end if
                end do
                !$omp end parallel do

                call vel_bc_staggered_CN2(u_star, v_star, w_star, bz_u_1, bz_u_nz, bz_v_1, bz_v_nz, bz_w_1, bz_w_nz, bc_z, 3)
                !$omp parallel do private(temp21) if (nx<40)
                do j=1,ny+2
                    temp21=transpose(u_star(:,j,:))
                    call dttrsb( F_low, F_d, F_up, temp21 )
                    u_star(:,j,:)=transpose(temp21)

                    if (j<=ny+1) then
                        temp21=transpose(v_star(:,j,:))
                        call dttrsb( F_low, F_d, F_up, temp21 )
                        v_star(:,j,:)=transpose(temp21)
                    end if

                    temp21=transpose(w_star(:,j,:))
                    call dttrsb( E_low, E_d, E_up, temp21 )
                    w_star(:,j,:)=transpose(temp21)
                end do
                !$omp end parallel do
            end if
            !call print_mat(u_star(:,:,2)-u_star_sub(:,:,2))

            !if (noise>0) then
            !    u_star = whiteNoise3(u_star, noise)
            !    v_star = whiteNoise3(v_star, noise)
            !    w_star = whiteNoise3(w_star, noise)
            !end if

            !RHS_poisson(2:ubound(RHS_poisson,1)-1,2:ubound(RHS_poisson,2)-1,2:ubound(RHS_poisson,3)-1) = &
            !    (diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx + &
            !    diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy + &
            !    diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz) / (dt0);
            RHS_poisson_internal = diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx
            RHS_poisson_internal = RHS_poisson_internal + diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy
            RHS_poisson_internal = RHS_poisson_internal + diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz
            RHS_poisson_internal = RHS_poisson_internal/dt0

            print *, [maxval(abs(u_star-u_star_sub)), maxval(abs(v_star-v_star_sub)), maxval(abs(w_star-w_star_sub)), maxval(abs(RHS_poisson_internal-RHS_poisson_sub(2:nxp-1,2:nyp-1,2:nzp-1)))]
            print *,'*****************************'

            temp31(2:nxp-1,2:nyp-1,2:nzp-1)=RHS_poisson_internal

            call pr_bc_staggered_modify_rhs(RHS_poisson_internal, bx_p_1(2:nyp-1,2:nzp-1), bx_p_nx(2:nyp-1,2:nzp-1), &
                by_p_1(2:nxp-1,2:nzp-1), by_p_ny(2:nxp-1,2:nzp-1), bz_p_1(2:nxp-1,2:nyp-1), bz_p_nz(2:nxp-1,2:nyp-1), &
                pbc_x, pbc_y, pbc_z, dx, dy, dz, dx2, dy2, dz2)
            print *,sum(RHS_poisson_internal)

            !!!!!!!!!! LU-based FD poisson solver !!!!!!!!!!
            if (LU_poisson) then
                RHS_poisson0=[RHS_poisson_internal]

                !DO i = 1, 64
                !    pt(i)%DUMMY = 0
                !END DO
                phase=33

                CALL SYSTEM_CLOCK(c1)
                call pardiso (pt, maxfct, mnum, mtype, phase, n, LHS_poisson%value, LHS_poisson%ia, LHS_poisson%ja, &
                    idum, nrhs, iparm, msglvl, RHS_poisson0, dp_vec, error)
                CALL SYSTEM_CLOCK(c2)
                print '("    Solve Poisson (LU decomp) completed: ", F8.4, " second")', (c2-c1)/system_clock_rate
                !print *, "**************************************"
                dp_lu(2:nxp-1,2:nyp-1,2:nzp-1)=reshape(dp_vec,([nx,ny,nz]))
                call pr_bc_staggered(dp_lu, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz, pbc_x, pbc_y, pbc_z, dx, dy, dz)
            end if

            !!!!!!!!!! FFT-based FD poisson solver !!!!!!!!!!
            if (FFT_poisson) then
                CALL SYSTEM_CLOCK(c1)
                !print *, "**************************************"
                !print *, "   Solve Poisson (FFT-based FD) start..."
                !dft_out_c=dcmplx(RHS_poisson_internal)

                !print *, "   Forward FFT..."
                status = DftiComputeForward(hand_f, RHS_poisson_internal(:,1,1), dft_out_c1(:,1,1))
                !status = DftiComputeForward(hand_f, dft_out_c(:,1,1))

                dft_out_c1=dft_out_c1/poisson_eigv
                !dft_out_c1(1,1,1)=0
                !print *, "   Backward FFT..."
                status = DftiComputeBackward(hand_b, dft_out_c1(:,1,1), dft_out_r1(:,1,1))
                !status = DftiComputeBackward(hand_b, dft_out_c(:,1,1))
                !dft_out_r=dble(dft_out_c)
                CALL SYSTEM_CLOCK(c2)
                print '("    Solve Poisson (FFT-based FD) completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

                dp(2:nxp-1,2:nyp-1,2:nzp-1)=dft_out_r1
                call pr_bc_staggered(dp, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz, pbc_x, pbc_y, pbc_z, dx, dy, dz)
                !print *, "**************************************"
            end if

            !!!!!!!!!! DCT-based FD poisson solver (pbc=1 and 4 are not implemented and tested) !!!!!!!!!!
            if (DCT_poisson) then
                CALL SYSTEM_CLOCK(c1)

                rhs_tt(1:nx,1:ny,1:nz)=RHS_poisson_internal

                call DCT_poisson_solver(rhs_tt, eig_tt, handle_x, handle_y, handle_z, &
                    ipar_x, ipar_y, ipar_z, dpar_x, dpar_y, dpar_z, nx, ny, nz, pbc_x, pbc_y, pbc_z, mean(dp_sub(2:nxp-1,2:nyp-1,2:nzp-1)))

                CALL SYSTEM_CLOCK(c2)
                print '("    Solve Poisson (DCT-based FD) completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

                dp(2:nxp-1,2:nyp-1,2:nzp-1)=rhs_tt(1:nx,1:ny,1:nz)
                call pr_bc_staggered(dp, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz, pbc_x, pbc_y, pbc_z, dx, dy, dz)
            end if
            !print *, mean(dp(2:nxp-1,2:nyp-1,2:nzp-1)), mean(dp_sub(2:nxp-1,2:nyp-1,2:nzp-1))
            !temp31=dp-dp_lu

            !OPEN(10, file="poisson_eq.dat", form="unformatted")
            !WRITE(10) [RHS_poisson]
            !WRITE(10) [dp]
            !WRITE(10) [dp_lu]
            !CLOSE(10)
            !dp=LHS_poisson\RHS_poisson(:);
            !dp=reshape(dp,nxp,nyp,nzp);

            !dp=dp_lu;
            !if (pbc_x==3 .and. pbc_y==3 .and. pbc_z==3) dp=dp-dp(1,1,1)+dp_sub(1,1,1)
            !$omp parallel sections
            !$omp section
            u=-dt0*diff(dp,1,1)/dx+u_star
            !$omp section
            v=-dt0*diff(dp,1,2)/dy+v_star
            !$omp section
            w=-dt0*diff(dp,1,3)/dz+w_star
            !$omp section
            p=dp_flag*p+dp
            !$omp end parallel sections

        end if
        !print *, p(3,4,5)
        CALL SYSTEM_CLOCK(c02)

        !!!!!!!!!!!!!!!!!!!!!!!!!! post-processing !!!!!!!!!!!!!!!!!!!!!!!!!!
        div=diff_old(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,1)/dx + &
            diff_old(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),1,2)/dy + &
            diff_old(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),1,3)/dz;

        rms_vel(1,t_step)=rms(u_sub);                             rms_vel(2,t_step)=rms(v_sub);                             rms_vel(3,t_step)=rms(w_sub);                             rms_vel(4,t_step)=rms(p_sub)
        err_vel(1,t_step)=maxval(abs(u-u_sub))/rms_vel(1,t_step); err_vel(2,t_step)=maxval(abs(v-v_sub))/rms_vel(2,t_step); err_vel(3,t_step)=maxval(abs(w-w_sub))/rms_vel(3,t_step); err_vel(4,t_step)=maxval(abs(p-p_sub))/rms_vel(4,t_step);
        err_vel(5,t_step)=rms(abs(u-u_sub))/rms_vel(1,t_step);    err_vel(6,t_step)=rms(abs(v-v_sub))/rms_vel(2,t_step);    err_vel(7,t_step)=rms(abs(w-w_sub))/rms_vel(3,t_step);    err_vel(8,t_step)=rms(abs(p-p_sub))/rms_vel(4,t_step);

        write(*,'("   MAX vel/pr error: ", 100g15.5)') err_vel(1:4,t_step)
        print '("Complete: ", F8.4, " second. MAX Div: ", E13.6)', (c02-c01)/system_clock_rate, maxval(abs(div))
        print *, "**************************************"

        if (save_output) then

            if (t_step==0) then
                write (string_var,'("err_file/err_vel_", A, ".h5")') trim(string_var2)

                call h5fcreate_f(string_var, H5F_ACC_EXCL_F, h5f_slice, status)

                !call h5ltset_attribute_int_f(h5f_slice, "/", "nx", [nx0], 1, status)
                !call h5ltset_attribute_double_f(h5f_slice, "/", "dt", [dt0], 1, status)
                !call h5ltset_attribute_double_f(h5f_slice, "/", "nu", [nu], 1, status)
                !call h5ltset_attribute_string_f(h5f_slice, "/", "time_scheme", trim(timescheme), status)
                !call h5ltset_attribute_int_f(h5f_slice, "/", "time_step", [t_step], 1, status)
                !call h5ltset_attribute_double_f(h5f_slice, "/", "time", [tGet], 1, status)
                !call h5ltset_attribute_int_f(h5f_slice, "/", "x0", [x0], 1, status)
                !call h5ltset_attribute_int_f(h5f_slice, "/", "nx0", [nx], 1, status)

                call h5ltset_attribute_string_f(h5f_slice, "/", "big_DNS_file", trim(big_DNS_file), status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "pbc_x", [pbc_x], 1, status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "using_Ustar", [bool2int(using_Ustar)], 1, status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "TOffset", [bool2int(TOffset)], 1, status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "restart", [bool2int(restart)], 1, status)
                call h5ltset_attribute_double_f(h5f_slice, "/", "noise", [noise], 1, status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "sub_tstep", [sub_tstep], 1, status)
            end if

            write (string_var,'("t_", F0.4)') tGet
            call h5gcreate_f(h5f_slice, string_var, h5g_sub, status)

            call h5ltmake_dataset_double_f(h5g_sub, 'u_sub',           3, [nx+1_8,ny+2_8,nz+2_8], u,           status)
            call h5ltmake_dataset_double_f(h5g_sub, 'v_sub',           3, [nx+2_8,ny+1_8,nz+2_8], v,           status)
            call h5ltmake_dataset_double_f(h5g_sub, 'w_sub',           3, [nx+2_8,ny+2_8,nz+1_8], w,           status)
            call h5ltmake_dataset_double_f(h5g_sub, 'p_sub',           3, [nx+2_8,ny+2_8,nz+2_8], p,           status)
            call h5ltmake_dataset_double_f(h5g_sub, 'u_star_sub',      3, [nx+1_8,ny+2_8,nz+2_8], u_star,      status)
            call h5ltmake_dataset_double_f(h5g_sub, 'v_star_sub',      3, [nx+2_8,ny+1_8,nz+2_8], v_star,      status)
            call h5ltmake_dataset_double_f(h5g_sub, 'w_star_sub',      3, [nx+2_8,ny+2_8,nz+1_8], w_star,      status)
            call h5ltmake_dataset_double_f(h5g_sub, 'dp_sub',          3, [nx+2_8,ny+2_8,nz+2_8], dp,          status)
            call h5ltmake_dataset_double_f(h5g_sub, 'RHS_poisson_sub', 3, [nx+2_8,ny+2_8,nz+2_8], temp31, status)

            call h5gclose_f( h5g_sub, status)

            if (tGet==t_end) then
                call h5fclose_f( h5f_slice, status)
            end if
        end if

    end do

    if (FFT_poisson) then
        status = DftiFreeDescriptor(hand_f)
        status = DftiFreeDescriptor(hand_b)
    end if

    if (DCT_poisson) then
        call free_trig_transform(handle_x,ipar_x,status)
        call free_trig_transform(handle_y,ipar_y,status)
        call free_trig_transform(handle_z,ipar_z,status)
    end if

    if (BC_interp) call h5fclose_f(h5f_bc_interp, status)
    call h5fclose_f(h5f_sub, status)
    call h5close_f(status)

    write (string_var,'("err_file/err_vel_", A, ".txt")') trim(string_var2)
    OPEN(10, file=string_var, form="formatted")
    write(10,'("time_step, rel_err_u, rel_err_v, rel_err_w, rel_err_p, rel_err_u_m, rel_err_v_m, rel_err_w_m, rel_err_p_m")')
    do j=0,time_length
        write(10,'(I0, ",", ES12.5E2, ",", ES12.5E2, ",", ES12.5E2, ",", ES12.5E2",", ES12.5E2",", ES12.5E2",", ES12.5E2",", ES12.5E2)') j, ( err_vel(i,j), i=1,8 )
    end do
    CLOSE(10)

    write (string_var,'("err_file/rms_vel_", A, "_", ES5.0E1, ".txt")') trim(timescheme), dt0
    OPEN(10, file=string_var, form="formatted")
    write(10,'("time_step, rms_u, rms_v, rms_w, rms_p")')
    do j=0,time_length
        write(10,'(I0, ",", ES12.5E2, ",", ES12.5E2, ",", ES12.5E2, ",", ES12.5E2)') j, ( rms_vel(i,j), i=1,4 )
    end do
    CLOSE(10)

    end subroutine re_simulation

