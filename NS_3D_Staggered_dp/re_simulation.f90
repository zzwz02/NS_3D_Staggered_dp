    subroutine re_simulation

    implicit none
    include 'mkl_pardiso.fi'
    !include 'fftw/fftw3.f'

    ! Variables
    real(8), parameter :: pi = 3.1415926535897932_8
    integer :: i=0,j=0,k=0,ll=0,mm=0

    logical, parameter :: init=.true., plot_output=.false.
    real(8), parameter :: Re=1.0d0, nu=0.002d0, t_end=20.0d0, dt0=0.004d0
    real(8) :: t_start
    !integer, parameter :: nx_file=256
    integer, parameter :: nx0=256, ny0=nx0, nz0=nx0, nxp0=nx0+2, nyp0=ny0+2, nzp0=nz0+2, sub_tstep=1
    real(8), parameter :: lx=2.0d0*pi, ly=2.0d0*pi, lz=2.0d0*pi, dx=lx/nx0, dy=ly/ny0, dz=lz/nz0, dx2=dx*dx, dy2=dy*dy, dz2=dz*dz
    real(8), parameter :: xu(nx0+1)=[(i*dx, i=0, nx0)],          yu(ny0+2)=[((i+0.5)*dy, i=-1, ny0)],   zu(nz0+2)=[((i+0.5)*dz, i=-1, nz0)]
    real(8), parameter :: xv(nx0+2)=[((i+0.5d0)*dx, i=-1, nx0)], yv(ny0+1)=[(i*dy, i=0, ny0)],          zv(nz0+2)=[((i+0.5d0)*dz, i=-1, nz0)]
    real(8), parameter :: xw(nx0+2)=[((i+0.5d0)*dx, i=-1, nx0)], yw(ny0+2)=[((i+0.5d0)*dy, i=-1, ny0)], zw(nz0+1)=[(i*dz, i=0, nz0)]
    real(8), parameter :: xp(nx0+2)=[((i+0.5d0)*dx, i=-1, nx0)], yp(ny0+2)=[((i+0.5d0)*dy, i=-1, ny0)], zp(nz0+2)=[((i+0.5d0)*dz, i=-1, nz0)]

    integer, parameter :: nx=32, ny=nx, nz=nx, nxp=nx+2, nyp=ny+2, nzp=nz+2
    integer, parameter :: x0=16, y0=16, z0=16
    integer, parameter :: idx_xu(nx+1)=[(x0+i, i=0, nx)],   idx_yu(ny+2)=[(y0+i, i=0, ny+1)], idx_zu(nz+2)=[(z0+i, i=0, nz+1)]
    integer, parameter :: idx_xv(nx+2)=[(x0+i, i=0, nx+1)], idx_yv(ny+1)=[(y0+i, i=0, ny)],   idx_zv(nz+2)=[(z0+i, i=0, nz+1)]
    integer, parameter :: idx_xw(nx+2)=[(x0+i, i=0, nx+1)], idx_yw(ny+2)=[(y0+i, i=0, ny+1)], idx_zw(nz+1)=[(z0+i, i=0, nz)]
    integer, parameter :: idx_xp(nx+2)=[(x0+i, i=0, nx+1)], idx_yp(ny+2)=[(y0+i, i=0, ny+1)], idx_zp(nz+2)=[(z0+i, i=0, nz+1)]

    logical, parameter :: LU_poisson=(nxp*nyp*nzp<=34**3)

    integer :: time_length
    real(8), dimension (:), allocatable :: time_array
    real(8) :: u(nx+1,ny+2,nz+2)=0,              v(nx+2,ny+1,nz+2)=0,              w(nx+2,ny+2,nz+1)=0,       p(nxp,nyp,nzp)=0
    real(8) :: u_star(nx+1,ny+2,nz+2)=0,         v_star(nx+2,ny+1,nz+2)=0,         w_star(nx+2,ny+2,nz+1)=0,  dp(nxp,nyp,nzp)=0
    real(8) :: rhs_x(nx+1,ny+2,nz+2)=0,          rhs_y(nx+2,ny+1,nz+2)=0,          rhs_z(nx+2,ny+2,nz+1)=0
    real(8) :: rhs_x_previous(nx+1,ny+2,nz+2)=0, rhs_y_previous(nx+2,ny+1,nz+2)=0, rhs_z_previous(nx+2,ny+2,nz+1)=0
    real(8) :: dpdx(nx+1,ny+2,nz+2)=0,           dpdy(nx+2,ny+1,nz+2)=0,           dpdz(nx+2,ny+2,nz+1)=0
    real(8) :: f_term_x=0,                       f_term_y=0,                       f_term_z=0
    !real(8) :: rhs_x_previous0(nx+1,ny+2,nz+2)=0,rhs_y_previous0(nx+2,ny+1,nz+2)=0,rhs_z_previous0(nx+2,ny+2,nz+1)=0
    !real(8) :: f_term_x(nx+1,ny+2,nz+2)=0,       f_term_y(nx+2,ny+1,nz+2)=0,       f_term_z(nx+2,ny+2,nz+1)=0

    real(8) :: RHS_poisson(nxp,nyp,nzp)=0, dp_lu(nxp,nyp,nzp)=0, RHS_poisson0(nxp*nyp*nzp)=0, dp_vec(nxp*nyp*nzp)=0
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
    real(8) :: RHS_poisson_sub(nx+2,ny+2,nz+2)=0

    !!!!!!!!!!!!!!! boundary conditions !!!!!!!!!!!!!!!
    real(8) :: bx_u_1(ny+2,nz+2)=0, bx_u_nx(ny+2,nz+2)=0, by_u_1(nx+1,nz+2)=0, by_u_ny(nx+1,nz+2)=0, bz_u_1(nx+1,ny+2)=0, bz_u_nz(nx+1,ny+2)=0
    real(8) :: bx_v_1(ny+1,nz+2)=0, bx_v_nx(ny+1,nz+2)=0, by_v_1(nx+2,nz+2)=0, by_v_ny(nx+2,nz+2)=0, bz_v_1(nx+2,ny+1)=0, bz_v_nz(nx+2,ny+1)=0
    real(8) :: bx_w_1(ny+2,nz+1)=0, bx_w_nx(ny+2,nz+1)=0, by_w_1(nx+2,nz+1)=0, by_w_ny(nx+2,nz+1)=0, bz_w_1(nx+2,ny+2)=0, bz_w_nz(nx+2,ny+2)=0
    real(8) :: bx_p_1(nyp,nzp)=0,   bx_p_nx(nyp,nzp)=0,   by_p_1(nxp,nzp)=0,   by_p_ny(nxp,nzp)=0,   bz_p_1(nxp,nyp)=0,   bz_p_nz(nxp,nyp)=0

    integer :: t_step, plot_step=20, slice=nz/2
    real(8) :: tGet
    integer, dimension (:), allocatable :: x_range0, y_range0, z_range0
    integer :: x_range1(nx)=[(i, i=2, nx+1)], y_range1(ny)=[(i, i=2, ny+1)], z_range1(nz)=[(i, i=2, nz+1)]
    type(csr), allocatable :: LHS_poisson
    !type(coo), allocatable :: LHS_poisson_coo
    !real(8), dimension (:,:), allocatable :: LHS_poisson_den

    !!!!!!!!!!!!!!! temperary variables !!!!!!!!!!!!!!!
    integer(8) :: sizeof_record, sizeof_record_sub, tempi1, tempi2
    !real(8) :: temp01, temp02, temp03!, temp04, temp05, temp06
    real(8), dimension (:), allocatable :: temp11!, temp12, temp13, temp14, temp15, temp16
    real(8), dimension (:,:), allocatable :: temp21!, temp22, temp23, temp24, temp25, temp26
    real(8), dimension (:,:,:), allocatable :: temp31, temp32, temp33!, temp34, temp35, temp36

    !!!!!!!!!!!!!!! INTEL mkl_pardiso !!!!!!!!!!!!!!!
    type(MKL_PARDISO_HANDLE) pt(64)
    integer :: maxfct=1, mnum=1, mtype=11, phase=13, n=nxp*nyp*nzp, idum(nxp*nyp*nzp), nrhs=1, iparm(64)=0, msglvl=0, error=0

    !!!!!!!!!!!!!!! INTEL mkl_dft !!!!!!!!!!!!!!!
    integer :: cstrides(4)=0, rstrides(4)=0
    type(DFTI_DESCRIPTOR), POINTER :: hand_f, hand_b
    !real(8) :: dft_out_r(nx,ny,nz), poisson_eigv(nx,ny,nz)=0
    !complex(8) :: dft_in_c(nx,ny,nz), dft_out_c(nx,ny,nz)
    real(8) :: dft_out_r1(nx,ny,nz), poisson_eigv(INT(nx/2.0)+1,ny,nz)=0
    complex(8) :: dft_out_c1(INT(nx/2.0)+1,ny,nz)
    integer :: status

    !!!!!!!!!!!!!!! post-processing !!!!!!!!!!!!!!!
    type(gpf):: gp
    character(len=1000) :: string_var

    !!!!!!!!!!!!!!! system_clock !!!!!!!!!!!!!!!
    REAL(8) :: system_clock_rate
    INTEGER :: c01,c02,c1,c2,cr,cm

    !!!!!!!!!!!!!!! re-simulation parameters !!!!!!!!!!!!!!!
    character(*), parameter :: timescheme="AB2-CN"
    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall); pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    integer, parameter :: bc_x=2, bc_y=bc_x, bc_z=bc_x, pbc_x=2, pbc_y=pbc_x, pbc_z=pbc_x
    logical, parameter :: using_Ustar=.true., TOffset=.false., restart=.false.
    real(8), parameter :: noise=0;
    real(8), dimension (:,:), allocatable :: err_vel, err_grad, rms_vel, rms_grad

    call OMP_set_dynamic(.true.)
    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)

    if (.not. TOffset) then
        t_start=0.0d0
    else
        t_start=1.0d0
    end if
    time_length=nint((t_end-t_start)/(dt0))
    allocate( time_array(0:time_length), err_vel(8,0:time_length), err_grad(24,0:time_length), rms_vel(4,0:time_length), rms_grad(12,0:time_length) )
    time_array=[(t_start+dt0*i, i=0, time_length)]

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
        allocate(A_low(nx), A_d(nx+1), A_up(nx), A_up2(nx-1), A_ipiv(n), B_low(nx+1), B_d(nx+2), B_up(nx+1), B_up2(nx), B_ipiv(n))
        allocate(C_low(nx), C_d(nx+1), C_up(nx), C_up2(nx-1), C_ipiv(n), D_low(nx+1), D_d(nx+2), D_up(nx+1), D_up2(nx), D_ipiv(n))
        allocate(E_low(nx), E_d(nx+1), E_up(nx), E_up2(nx-1), E_ipiv(n), F_low(nx+1), F_d(nx+2), F_up(nx+1), F_up2(nx), F_ipiv(n))

        call mat_CN_tridiagonal(nx, bc_x, dx2, dx, A_low, A_d, A_up, B_low, B_d, B_up, dt0, nu)
        call mat_CN_tridiagonal(ny, bc_y, dy2, dy, C_low, C_d, C_up, D_low, D_d, D_up, dt0, nu)
        call mat_CN_tridiagonal(nz, bc_z, dz2, dz, E_low, E_d, E_up, F_low, F_d, F_up, dt0, nu)

        call gttrf( A_low, A_d, A_up, A_up2, A_ipiv )
        call gttrf( B_low, B_d, B_up, B_up2, B_ipiv )
        call gttrf( C_low, C_d, C_up, C_up2, C_ipiv )
        call gttrf( D_low, D_d, D_up, D_up2, D_ipiv )
        call gttrf( E_low, E_d, E_up, E_up2, E_ipiv )
        call gttrf( F_low, F_d, F_up, F_up2, F_ipiv )
    end if

    if (LU_poisson) then
        LHS_poisson=Poisson_LHS_staggered(nxp, nyp, nzp, dx2, dy2, dz2, pbc_x, pbc_y, pbc_z, dx, dy, dz)
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
        !iparm(2) = 2 ! fill-in reordering from METIS
        !iparm(5) = 2 ! no user fill-in reducing permutation
        !iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        !iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
        !iparm(19) = -1 ! Output: Mflops for LU factorization

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
    poisson_eigv(1,1,1)=1.0d0
    poisson_eigv=-poisson_eigv*4.0d0

    sizeof_record = (nx0+1)*(ny0+2)*(nz0+2) + (nx0+2)*(ny0+1)*(nz0+2) + (nx0+2)*(ny0+2)*(nz0+1) + (nx0+2)*(ny0+2)*(nz0+2) + &
        (nx0+1)*(ny0+2)*(nz0+2) + (nx0+2)*(ny0+1)*(nz0+2) + (nx0+2)*(ny0+2)*(nz0+1)
    sizeof_record_sub=size(u_sub)+size(v_sub)+size(w_sub)+size(p_sub)+size(u_star_sub)+size(v_star_sub)+size(w_star_sub)+size(dp_sub)+size(RHS_poisson_sub)
    if (timescheme=="AB2-CN") then
        sizeof_record_sub=sizeof_record_sub+ &
            size(bx_u_1)+size(bx_u_nx)+size(bx_v_1)+size(bx_v_nx)+size(bx_w_1)+size(bx_w_nx)+ &
            size(by_u_1)+size(by_u_ny)+size(by_v_1)+size(by_v_ny)+size(by_w_1)+size(by_w_ny)+ &
            size(bz_u_1)+size(bz_u_nz)+size(bz_v_1)+size(bz_v_nz)+size(bz_w_1)+size(bz_w_nz) !+size(u_sub)+size(v_sub)+size(w_sub)
    end if

    do t_step=0,time_length
        tGet=time_array(t_step)
        print *,''
        print '(" t_step ", I6, "        tGet ", F7.4)', t_step, tGet
        CALL SYSTEM_CLOCK(c01)

        if (t_step==0) then
            if (allocated(temp11)) deallocate(temp11)
            allocate(temp11(sizeof_record))
            write (string_var,'(A, "_result\HIT_", I0,"^3_decay_", ES5.0E1, "_", A ,"_dp_t_" , F0.4, ".dat")') trim(timescheme), nx0, dt0, trim(timescheme), tGet
            open(20, file=string_var, form='unformatted', status='old', action='read', &
                access='direct', recl=sizeof_record*2) !variables are double precision
            read(20, rec=1) temp11
            close(20)

            tempi1=1;        tempi2=(nx0+1)*(ny0+2)*(nz0+2);        temp31=reshape(temp11(tempi1:tempi2), [nx0+1,ny0+2,nz0+2]); u_sub=temp31(idx_xu,idx_yu,idx_zu)
            tempi1=tempi2+1; tempi2=tempi2+(nx0+2)*(ny0+1)*(nz0+2); temp31=reshape(temp11(tempi1:tempi2), [nx0+2,ny0+1,nz0+2]); v_sub=temp31(idx_xv,idx_yv,idx_zv)
            tempi1=tempi2+1; tempi2=tempi2+(nx0+2)*(ny0+2)*(nz0+1); temp31=reshape(temp11(tempi1:tempi2), [nx0+2,ny0+2,nz0+1]); w_sub=temp31(idx_xw,idx_yw,idx_zw)
            tempi1=tempi2+1; tempi2=tempi2+(nx0+2)*(ny0+2)*(nz0+2); temp31=reshape(temp11(tempi1:tempi2), [nx0+2,ny0+2,nz0+2]); p_sub=temp31(idx_xp,idx_yp,idx_zp)
            tempi1=tempi2+1; tempi2=tempi2+(nx0+1)*(ny0+2)*(nz0+2); temp31=reshape(temp11(tempi1:tempi2), [nx0+1,ny0+2,nz0+2]); rhs_x_previous=temp31(idx_xu,idx_yu,idx_zu)
            tempi1=tempi2+1; tempi2=tempi2+(nx0+2)*(ny0+1)*(nz0+2); temp31=reshape(temp11(tempi1:tempi2), [nx0+2,ny0+1,nz0+2]); rhs_y_previous=temp31(idx_xv,idx_yv,idx_zv)
            tempi1=tempi2+1; tempi2=tempi2+(nx0+2)*(ny0+2)*(nz0+1); temp31=reshape(temp11(tempi1:tempi2), [nx0+2,ny0+2,nz0+1]); rhs_z_previous=temp31(idx_xw,idx_yw,idx_zw)
            deallocate(temp11, temp31)
            u=u_sub; v=v_sub; w=w_sub; p=p_sub;
            !call TGV(xu, yu, zu, 0.0d0, nu, u)
            !call TGV(xv, yv, zv, 0.0d0, nu, v=v)
            !call TGV(xp, yp, zp, 0.0d0, nu, p=p)
        else
            if (allocated(temp11) .and. size(temp11)/=sizeof_record_sub) then
                deallocate(temp11)
            end if
            if (.not. allocated(temp11)) allocate(temp11(sizeof_record_sub))

            write (string_var,'(A, "_result\HIT_", I0,"^3_decay_", ES5.0E1, "_", A , "_dp_x0_", I0, "_nx0_", I0, "_t_", F0.4, ".dat")') trim(timescheme), nx0, dt0, trim(timescheme), x0, nx, tGet
            open(20, file=string_var, form='unformatted', status='old', action='read', &
                access='direct', recl=sizeof_record_sub*2) !variables are double precision
            read(20, rec=1) temp11
            close(20)

            tempi1=1;        tempi2=size(u_sub);             u_sub=reshape(temp11(tempi1:tempi2), [nx+1,ny+2,nz+2]);
            tempi1=tempi2+1; tempi2=tempi2+size(v_sub);      v_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+1,nz+2]);
            tempi1=tempi2+1; tempi2=tempi2+size(w_sub);      w_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+1]);
            tempi1=tempi2+1; tempi2=tempi2+size(p_sub);      p_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+2]);
            tempi1=tempi2+1; tempi2=tempi2+size(u_star_sub); u_star_sub=reshape(temp11(tempi1:tempi2), [nx+1,ny+2,nz+2]);
            tempi1=tempi2+1; tempi2=tempi2+size(v_star_sub); v_star_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+1,nz+2]);
            tempi1=tempi2+1; tempi2=tempi2+size(w_star_sub); w_star_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+1]);
            tempi1=tempi2+1; tempi2=tempi2+size(dp_sub);     dp_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+2]);
            tempi1=tempi2+1; tempi2=tempi2+size(RHS_poisson_sub);     RHS_poisson_sub=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+2]);
            if (timescheme=="AB2-CN") then
                tempi1=tempi2+1; tempi2=tempi2+size(bx_u_1);  bx_u_1 =reshape(temp11(tempi1:tempi2), [ny+2,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bx_u_nx); bx_u_nx=reshape(temp11(tempi1:tempi2), [ny+2,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bx_v_1);  bx_v_1 =reshape(temp11(tempi1:tempi2), [ny+1,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bx_v_nx); bx_v_nx=reshape(temp11(tempi1:tempi2), [ny+1,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bx_w_1);  bx_w_1 =reshape(temp11(tempi1:tempi2), [ny+2,nz+1]);
                tempi1=tempi2+1; tempi2=tempi2+size(bx_w_nx); bx_w_nx=reshape(temp11(tempi1:tempi2), [ny+2,nz+1]);

                tempi1=tempi2+1; tempi2=tempi2+size(by_u_1);  by_u_1 =reshape(temp11(tempi1:tempi2), [nx+1,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(by_u_ny); by_u_ny=reshape(temp11(tempi1:tempi2), [nx+1,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(by_v_1);  by_v_1 =reshape(temp11(tempi1:tempi2), [nx+2,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(by_v_ny); by_v_ny=reshape(temp11(tempi1:tempi2), [nx+2,nz+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(by_w_1);  by_w_1 =reshape(temp11(tempi1:tempi2), [nx+2,nz+1]);
                tempi1=tempi2+1; tempi2=tempi2+size(by_w_ny); by_w_ny=reshape(temp11(tempi1:tempi2), [nx+2,nz+1]);

                tempi1=tempi2+1; tempi2=tempi2+size(bz_u_1);  bz_u_1 =reshape(temp11(tempi1:tempi2), [nx+1,ny+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bz_u_nz); bz_u_nz=reshape(temp11(tempi1:tempi2), [nx+1,ny+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bz_v_1);  bz_v_1 =reshape(temp11(tempi1:tempi2), [nx+2,ny+1]);
                tempi1=tempi2+1; tempi2=tempi2+size(bz_v_nz); bz_v_nz=reshape(temp11(tempi1:tempi2), [nx+2,ny+1]);
                tempi1=tempi2+1; tempi2=tempi2+size(bz_w_1);  bz_w_1 =reshape(temp11(tempi1:tempi2), [nx+2,ny+2]);
                tempi1=tempi2+1; tempi2=tempi2+size(bz_w_nz); bz_w_nz=reshape(temp11(tempi1:tempi2), [nx+2,ny+2]);
                
                !tempi1=tempi2+1; tempi2=tempi2+size(u_sub);   temp31=reshape(temp11(tempi1:tempi2), [nx+1,ny+2,nz+2]);
                !tempi1=tempi2+1; tempi2=tempi2+size(v_sub);   temp32=reshape(temp11(tempi1:tempi2), [nx+2,ny+1,nz+2]);
                !tempi1=tempi2+1; tempi2=tempi2+size(w_sub);   temp33=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+1]);
            end if

            if (timescheme=="AB2-CN") then
                call get_pr_bc(dp_sub, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)
            else if (using_Ustar) then
                call get_velpr_bc(u_star_sub, v_star_sub, w_star_sub, dp_sub, bc_x, bc_y, bc_z, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, &
                    bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz, bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz, &
                    bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)
            else
                call get_velpr_bc(u_sub, v_sub, w_sub, dp_sub, bc_x, bc_y, bc_z, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, &
                    bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz, bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz, &
                    bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)
            end if

            !!! Convection
            call cal_conv(u, v, w, bc_x, bc_y, bc_z, dx, dy, dz, conv_x, conv_y, conv_z)

            !!! Diffusion
            call cal_diff(u, v, w, bc_x, bc_y, bc_z, dx2, dy2, dz2, diff_x, diff_y, diff_z)

            !!! Time-advancement
            f_term_x=0; f_term_y=0; f_term_z=0;

            dpdx=diff(p,1,1)/dx;
            dpdy=diff(p,1,2)/dy;
            dpdz=diff(p,1,3)/dz;

            if (timescheme=="Euler") then
                ! diff terms + conv terms
                rhs_x(x_range0,y_range1,z_range1)=nu*diff_x-conv_x;
                rhs_y(x_range1,y_range0,z_range1)=nu*diff_y-conv_y;
                rhs_z(x_range1,y_range1,z_range0)=nu*diff_z-conv_z;

                u_star=dt0*(1.0d0*rhs_x-1.0d0*dpdx+1.0d0*f_term_x)+u;
                v_star=dt0*(1.0d0*rhs_y-1.0d0*dpdy+1.0d0*f_term_y)+v;
                w_star=dt0*(1.0d0*rhs_z-1.0d0*dpdz+1.0d0*f_term_z)+w;

                rhs_x_previous=0;
                rhs_y_previous=0;
                rhs_z_previous=0;

                call vel_bc_staggered(u_star,v_star,w_star,&
                    bx_u_1,bx_u_nx,by_u_1,by_u_ny,bz_u_1,bz_u_nz,&
                    bx_v_1,bx_v_nx,by_v_1,by_v_ny,bz_v_1,bz_v_nz,&
                    bx_w_1,bx_w_nx,by_w_1,by_w_ny,bz_w_1,bz_w_nz,&
                    bc_x,bc_y,bc_z);
            else if (timescheme=="AB2") then
                ! diff terms + conv terms
                rhs_x(x_range0,y_range1,z_range1)=nu*diff_x-conv_x;
                rhs_y(x_range1,y_range0,z_range1)=nu*diff_y-conv_y;
                rhs_z(x_range1,y_range1,z_range0)=nu*diff_z-conv_z;

                ! prediction
                if ((.not. restart .and. t_step==1) .or. (.not. TOffset .and. t_step==1)) then
                    u_star=dt0*(1.0d0*rhs_x-1.0d0*dpdx+1.0d0*f_term_x)+u;
                    v_star=dt0*(1.0d0*rhs_y-1.0d0*dpdy+1.0d0*f_term_y)+v;
                    w_star=dt0*(1.0d0*rhs_z-1.0d0*dpdz+1.0d0*f_term_z)+w;
                else
                    u_star=dt0*(1.5d0*rhs_x-0.5d0*rhs_x_previous-1.0d0*dpdx+1.0d0*f_term_x)+u;
                    v_star=dt0*(1.5d0*rhs_y-0.5d0*rhs_y_previous-1.0d0*dpdy+1.0d0*f_term_y)+v;
                    w_star=dt0*(1.5d0*rhs_z-0.5d0*rhs_z_previous-1.0d0*dpdz+1.0d0*f_term_z)+w;

                end if

                rhs_x_previous=rhs_x;
                rhs_y_previous=rhs_y;
                rhs_z_previous=rhs_z;

                call vel_bc_staggered(u_star,v_star,w_star,&
                    bx_u_1,bx_u_nx,by_u_1,by_u_ny,bz_u_1,bz_u_nz,&
                    bx_v_1,bx_v_nx,by_v_1,by_v_ny,bz_v_1,bz_v_nz,&
                    bx_w_1,bx_w_nx,by_w_1,by_w_ny,bz_w_1,bz_w_nz,&
                    bc_x,bc_y,bc_z);
            else if (timescheme=="AB2-CN") then
                ! prediction
                !if (t_step==1 .and. init) then
                if ( (t_step==1 .or. t_step==nint( (dt0-t_start)/dt0) ) .and. init) then
                    ! diff terms + conv terms
                    rhs_x(x_range0,y_range1,z_range1)=0.5d0*nu*diff_x-conv_x;
                    rhs_y(x_range1,y_range0,z_range1)=0.5d0*nu*diff_y-conv_y;
                    rhs_z(x_range1,y_range1,z_range0)=0.5d0*nu*diff_z-conv_z;

                    u_star=dt0*(rhs_x-1.0d0*dpdx+1.0d0*f_term_x)+u;
                    v_star=dt0*(rhs_y-1.0d0*dpdy+1.0d0*f_term_y)+v;
                    w_star=dt0*(rhs_z-1.0d0*dpdz+1.0d0*f_term_z)+w;
                else
                    ! diff terms + conv terms
                    rhs_x(x_range0,y_range1,z_range1)=0.5d0*nu*diff_x-1.5d0*conv_x;
                    rhs_y(x_range1,y_range0,z_range1)=0.5d0*nu*diff_y-1.5d0*conv_y;
                    rhs_z(x_range1,y_range1,z_range0)=0.5d0*nu*diff_z-1.5d0*conv_z;

                    u_star=dt0*(rhs_x-0.5d0*rhs_x_previous-1.0d0*dpdx+1.0d0*f_term_x)+u;
                    v_star=dt0*(rhs_y-0.5d0*rhs_y_previous-1.0d0*dpdy+1.0d0*f_term_y)+v;
                    w_star=dt0*(rhs_z-0.5d0*rhs_z_previous-1.0d0*dpdz+1.0d0*f_term_z)+w;
                end if

                call vel_bc_staggered_CN2(u_star, v_star, w_star, bx_u_1, bx_u_nx, bx_v_1, bx_v_nx, bx_w_1, bx_w_nx, bc_x, 1)
                !!$omp parallel do
                do k=1,nz+2
                    call gttrs( A_low, A_d, A_up, A_up2, u_star(:,:,k), A_ipiv )
                    call gttrs( B_low, B_d, B_up, B_up2, v_star(:,:,k), B_ipiv )
                    if (k<=nz+1) call gttrs( B_low, B_d, B_up, B_up2, w_star(:,:,k), B_ipiv )
                end do
                !!$omp end parallel do
                
                call vel_bc_staggered_CN2(u_star, v_star, w_star, by_u_1, by_u_ny, by_v_1, by_v_ny, by_w_1, by_w_ny, bc_y, 2)
                !!$omp parallel do
                do k=1,nz+2
                    temp21=transpose(u_star(:,:,k))
                    call gttrs( D_low, D_d, D_up, D_up2, temp21, D_ipiv )
                    u_star(:,:,k)=transpose(temp21)

                    temp21=transpose(v_star(:,:,k))
                    call gttrs( C_low, C_d, C_up, C_up2, temp21, C_ipiv )
                    v_star(:,:,k)=transpose(temp21)

                    if (k<=nz+1) then
                        temp21=transpose(w_star(:,:,k))
                        call gttrs( D_low, D_d, D_up, D_up2, temp21, D_ipiv )
                        w_star(:,:,k)=transpose(temp21)
                    end if
                end do
                !!$omp end parallel do

                call vel_bc_staggered_CN2(u_star, v_star, w_star, bz_u_1, bz_u_nz, bz_v_1, bz_v_nz, bz_w_1, bz_w_nz, bc_z, 3)
                !!$omp parallel do
                do j=1,ny+2
                    temp21=transpose(u_star(:,j,:))
                    call gttrs( F_low, F_d, F_up, F_up2, temp21, F_ipiv )
                    u_star(:,j,:)=transpose(temp21)

                    if (j<=ny+1) then
                        temp21=transpose(v_star(:,j,:))
                        call gttrs( F_low, F_d, F_up, F_up2, temp21, F_ipiv )
                        v_star(:,j,:)=transpose(temp21)
                    end if

                    temp21=transpose(w_star(:,j,:))
                    call gttrs( E_low, E_d, E_up, E_up2, temp21, E_ipiv )
                    w_star(:,j,:)=transpose(temp21)
                end do
                !!$omp end parallel do

                rhs_x_previous(x_range0,y_range1,z_range1)=-conv_x;
                rhs_y_previous(x_range1,y_range0,z_range1)=-conv_y;
                rhs_z_previous(x_range1,y_range1,z_range0)=-conv_z;
            end if

            !call print_mat(u_star(:,:,2)-u_star_sub(:,:,2))

            !RHS_poisson(2:ubound(RHS_poisson,1)-1,2:ubound(RHS_poisson,2)-1,2:ubound(RHS_poisson,3)-1) = &
            !    (diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx + &
            !    diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy + &
            !    diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz) / (dt0);
            RHS_poisson_internal = diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx
            RHS_poisson_internal = RHS_poisson_internal + diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy
            RHS_poisson_internal = RHS_poisson_internal + diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz
            RHS_poisson_internal = RHS_poisson_internal/dt0
            RHS_poisson(2:nxp-1,2:nyp-1,2:nzp-1)=RHS_poisson_internal

            print *, [maxval(abs(u_star-u_star_sub)), maxval(abs(v_star-v_star_sub)), maxval(abs(w_star-w_star_sub)), maxval(abs(RHS_poisson_internal-RHS_poisson_sub(2:nxp-1,2:nyp-1,2:nzp-1)))]
            print *,'*****************************'

            if (LU_poisson) then
                !RHS_poisson(2:nxp-1,2:nyp-1,2:nzp-1)=RHS_poisson_internal
                if (pbc_x==1) then
                    RHS_poisson(1,:,:)=0;      RHS_poisson(nxp,:,:)=0;
                else if (pbc_x==2 .or. pbc_x==4) then
                    RHS_poisson(1,:,:)=bx_p_1; RHS_poisson(nxp,:,:)=bx_p_nx;
                else if (pbc_x==3) then
                    RHS_poisson(1,:,:)=bx_p_1; RHS_poisson(nxp,:,:)=bx_p_nx;
                end if
                if (pbc_y==1) then
                    RHS_poisson(:,1,:)=0;      RHS_poisson(:,nyp,:)=0;
                else if (pbc_y==2 .or. pbc_y==4) then
                    RHS_poisson(:,1,:)=by_p_1; RHS_poisson(:,nyp,:)=by_p_ny;
                else if (pbc_y==3) then
                    RHS_poisson(:,1,:)=by_p_1; RHS_poisson(:,nyp,:)=by_p_ny;
                end if
                if (pbc_z==1) then
                    RHS_poisson(:,:,1)=0;      RHS_poisson(:,:,nzp)=0;
                else if (pbc_y==2 .or. pbc_y==4) then
                    RHS_poisson(:,:,1)=bz_p_1; RHS_poisson(:,:,nzp)=bz_p_nz;
                else if (pbc_z==3) then
                    RHS_poisson(:,:,1)=bz_p_1; RHS_poisson(:,:,nzp)=bz_p_nz;
                end if

                RHS_poisson0=[RHS_poisson]

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
                dp_lu=reshape(dp_vec,([nxp,nyp,nzp]))
            end if

            CALL SYSTEM_CLOCK(c1)
            !print *, "**************************************"
            !print *, "   Solve Poisson (FFT-based FD) start..."
            !dft_out_c=dcmplx(RHS_poisson_internal)

            !print *, "   Forward FFT..."
            status = DftiComputeForward(hand_f, RHS_poisson_internal(:,1,1), dft_out_c1(:,1,1))
            !status = DftiComputeForward(hand_f, dft_out_c(:,1,1))

            dft_out_c1=dft_out_c1/poisson_eigv
            dft_out_c1(1,1,1)=0
            !print *, "   Backward FFT..."
            status = DftiComputeBackward(hand_b, dft_out_c1(:,1,1), dft_out_r1(:,1,1))
            !status = DftiComputeBackward(hand_b, dft_out_c(:,1,1))
            !dft_out_r=dble(dft_out_c)
            dp(2:nxp-1,2:nyp-1,2:nzp-1)=dft_out_r1
            dp(1,:,:)=dp(nxp-1,:,:); dp(nxp,:,:)=dp(2,:,:)
            dp(:,1,:)=dp(:,nyp-1,:); dp(:,nyp,:)=dp(:,2,:)
            dp(:,:,1)=dp(:,:,nzp-1); dp(:,:,nzp)=dp(:,:,2)

            CALL SYSTEM_CLOCK(c2)
            print '("    Solve Poisson (FFT-based FD) completed: ", F8.4, " second")', (c2-c1)/system_clock_rate
            !print *, "**************************************"

            !print *, maxval(abs(dp-dp(1,1,1)-dp_lu+dp_lu(1,1,1)))
            !temp31=dp-dp_lu

            !OPEN(10, file="poisson_eq.dat", form="unformatted")
            !WRITE(10) [RHS_poisson]
            !WRITE(10) [dp]
            !WRITE(10) [dp_lu]
            !CLOSE(10)
            !dp=LHS_poisson\RHS_poisson(:);
            !dp=reshape(dp,nxp,nyp,nzp);

            dp=dp_lu; !dp=dp-dp(1,1,1)+dp_sub(1,1,1)
            u=-dt0*diff(dp,1,1)/dx+u_star
            v=-dt0*diff(dp,1,2)/dy+v_star
            w=-dt0*diff(dp,1,3)/dz+w_star
            p=p+dp

        end if
        !print *, p(3,4,5)
        CALL SYSTEM_CLOCK(c02)

        !!!!!!!!!!!!!!!!!!!!!!!!!! post-processing !!!!!!!!!!!!!!!!!!!!!!!!!!
        div=diff_old(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,1)/dx + &
            diff_old(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),1,2)/dy + &
            diff_old(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),1,3)/dz;

        !real(8), dimension (:,:), allocatable :: err_vel, err_grad, rms_vel, rms_grad
        !allocate( time_array(0:time_length), err_vel(0:time_length,8), err_grad(0:time_length,24), rms_vel(0:time_length,4), rms_grad(0:time_length,12) )
        rms_vel(1,t_step)=rms(u_sub);                               rms_vel(2,t_step)=rms(v_sub);                               rms_vel(3,t_step)=rms(w_sub);                               rms_vel(4,t_step)=rms(p_sub)
        err_vel(1,t_step)=maxval(abs(u-u_sub))/rms_vel(1,t_step);   err_vel(2,t_step)=maxval(abs(v-v_sub))/rms_vel(2,t_step);   err_vel(3,t_step)=maxval(abs(w-w_sub))/rms_vel(3,t_step);   err_vel(4,t_step)=maxval(abs(dp-dp_sub))/rms_vel(4,t_step)
        err_vel(5,t_step)=mean(u-u_sub)/rms_vel(1,t_step);          err_vel(6,t_step)=mean(v-v_sub)/rms_vel(2,t_step);          err_vel(7,t_step)=mean(w-w_sub)/rms_vel(3,t_step);          err_vel(8,t_step)=mean(p-p_sub)/rms_vel(4,t_step)

        !rms_grad(1,t_step)=rms(u_sub);               rms_grad(2,t_step)=rms(v_sub);               rms_grad(3,t_step)=rms(w_sub);               rms_grad(4,t_step)=rms(p_sub)
        !err_grad(1,t_step)=maxval(u-u_sub);          err_grad(2,t_step)=maxval(v-v_sub);          err_grad(3,t_step)=maxval(w-w_sub);          err_grad(4,t_step)=maxval(p-p_sub)
        !err_grad(5,t_step)=mean(u-u_sub)/rms(u_sub); err_grad(6,t_step)=mean(v-v_sub)/rms(v_sub); err_grad(7,t_step)=mean(w-w_sub)/rms(w_sub); err_grad(8,t_step)=mean(p-p_sub)/rms(p_sub)

        write(*,'("   MAX vel/pr error: ", 100g15.5)') err_vel(1:4,t_step)
        print '("Complete: ", F8.4, " second. MAX Div: ", E13.6)', (c02-c01)/system_clock_rate, maxval(abs(div))
        print *, "**************************************"

        if (mod(t_step,plot_step)==0 .and. plot_output) then
            plot_block: block

                !real(8), dimension (:,:), allocatable :: xx, yy

                !call meshgrid( xx, yy, yu(2:ubound(yu,1)-1), xu )
                write (string_var,'("re-sim_t_",i0.4)') t_step
                call gp%output_filename(string_var)
                call gp%multiplot(2,2)

                write (string_var,'("U @ ",i7)') t_step
                call gp%title(string_var)
                !call gp%xlabel('x1, ...')
                !call gp%ylabel('y1, ...')
                call gp%contour(u(:,2:ubound(u,2)-1,slice), palette='jet')

                !write (string_var,'("V")') t_step
                call gp%title("V")
                !call gp%xlabel('x1, ...')
                !call gp%ylabel('y1, ...')
                call gp%contour(v(2:ubound(v,1)-1,:,slice), palette='jet')

                !write (string_var,'("W")') t_step
                call gp%title("W")
                !call gp%xlabel('x1, ...')
                !call gp%ylabel('y1, ...')
                call gp%contour(w(2:ubound(w,1)-1,2:ubound(w,2)-1,slice), palette='jet')

                !write (string_var,'("P")') t_step
                call gp%title("P")
                !call gp%xlabel('x1, ...')
                !call gp%ylabel('y1, ...')
                call gp%contour(p(2:ubound(p,1)-1,2:ubound(p,2)-1,slice), palette='jet')
                if (gp%hasfileopen) call draw_now(gp)

                !subplot(2,3,1)
                !temp=squeeze(u(:,2:end-1,slice));
                !contourf(temp',[-1e99 linspace(min(temp(:)),max(temp(:)),21)],'LineStyle','none'); colorbar;
                !caxis([min(temp(:)),max(temp(:))]);
                !title(['U @ ',num2str(t_step)]);
                !set(gca,'FontSize',18);
                !
                !subplot(2,3,2)
                !temp=squeeze(v(2:end-1,:,slice));
                !contourf(temp',[-1e99 linspace(min(temp(:)),max(temp(:)),21)],'LineStyle','none'); colorbar;
                !caxis([min(temp(:)),max(temp(:))]);
                !title(['V']);
                !set(gca,'FontSize',18);
                !
                !subplot(2,3,3)
                !temp=squeeze(w(2:end-1,2:end-1,slice));
                !contourf(temp',[-1e99 linspace(min(temp(:)),max(temp(:)),21)],'LineStyle','none'); colorbar;
                !caxis([min(temp(:)),max(temp(:))]);
                !title(['W']);
                !set(gca,'FontSize',18);
                !
                !subplot(2,3,4)
                !temp=squeeze(p(2:end-1,2:end-1,slice));
                !contourf(temp',[-1e99 linspace(min(temp(:)),max(temp(:)),21)],'LineStyle','none'); colorbar;
                !caxis([min(temp(:)) max(temp(:))]);
                !title(['P']);
                !set(gca,'FontSize',18);
                !
                !drawnow

            end block plot_block
        end if


    end do

    status = DftiFreeDescriptor(hand_f)
    status = DftiFreeDescriptor(hand_b)

    OPEN(10, file="err_vel.txt", form="formatted")
    write(10,'("rel_err")')
    do j=0,time_length
        write(10,'(100g15.5)') ( err_vel(i,j), i=1,4 )
    end do
    write(10,'("***************************************")')
    write(10,'("rms_val")')
    do j=0,time_length
        write(10,'(100g15.5)') ( rms_vel(i,j), i=1,4 )
    end do
    CLOSE(10)

    end subroutine re_simulation

