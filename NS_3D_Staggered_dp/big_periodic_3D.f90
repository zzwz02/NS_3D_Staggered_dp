    subroutine big_periodic_3D

    implicit none

    ! Variables
    real(8), parameter :: pi = 3.1415926535897932_8
    integer :: i=0,j=0,k=0,ll=0,mm=0

    logical, parameter :: init=.true.
    real(8), parameter :: Re=1.0d0, nu=0.002d0, t_end=10.0d0!, dt0=0.004d0 !!!dt<=0.004 to be stable
    real(8), parameter :: t_start=0.0d0
    !integer, parameter :: nx_file=256
    integer, parameter :: nx=256, ny=nx, nz=nx, nxp=nx+2, nyp=ny+2, nzp=nz+2, sub_tstep=1
    real(8), parameter :: lx=2.0d0*pi, ly=2.0d0*pi, lz=2.0d0*pi, dx=lx/nx, dy=ly/ny, dz=lz/nz, dx2=dx*dx, dy2=dy*dy, dz2=dz*dz
    real(8), parameter :: xu(nx+1)=[(i*dx, i=0, nx)],          yu(ny+2)=[((i+0.5)*dy, i=-1, ny)],   zu(nz+2)=[((i+0.5)*dz, i=-1, nz)]
    real(8), parameter :: xv(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yv(ny+1)=[(i*dy, i=0, ny)],          zv(nz+2)=[((i+0.5d0)*dz, i=-1, nz)]
    real(8), parameter :: xw(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yw(ny+2)=[((i+0.5d0)*dy, i=-1, ny)], zw(nz+1)=[(i*dz, i=0, nz)]
    real(8), parameter :: xp(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yp(ny+2)=[((i+0.5d0)*dy, i=-1, ny)], zp(nz+2)=[((i+0.5d0)*dz, i=-1, nz)]

    integer, parameter :: nx0=32, ny0=nx0, nz0=nx0
    integer, parameter :: x0=16, y0=16, z0=16
    integer, parameter :: idx_xu(nx0+1)=[(x0+i, i=0, nx0)],   idx_yu(ny0+2)=[(y0+i, i=0, ny0+1)], idx_zu(nz0+2)=[(z0+i, i=0, nz0+1)]
    integer, parameter :: idx_xv(nx0+2)=[(x0+i, i=0, nx0+1)], idx_yv(ny0+1)=[(y0+i, i=0, ny0)],   idx_zv(nz0+2)=[(z0+i, i=0, nz0+1)]
    integer, parameter :: idx_xw(nx0+2)=[(x0+i, i=0, nx0+1)], idx_yw(ny0+2)=[(y0+i, i=0, ny0+1)], idx_zw(nz0+1)=[(z0+i, i=0, nz0)]
    integer, parameter :: idx_xp(nx0+2)=[(x0+i, i=0, nx0+1)], idx_yp(ny0+2)=[(y0+i, i=0, ny0+1)], idx_zp(nz0+2)=[(z0+i, i=0, nz0+1)]

    real(8) :: u(nx+1,ny+2,nz+2)=0,              v(nx+2,ny+1,nz+2)=0,              w(nx+2,ny+2,nz+1)=0,       p(nxp,nyp,nzp)=0
    real(8) :: u_star(nx+1,ny+2,nz+2)=0,         v_star(nx+2,ny+1,nz+2)=0,         w_star(nx+2,ny+2,nz+1)=0,  dp(nxp,nyp,nzp)=0
    !real(8) :: rhs_x(nx+1,ny+2,nz+2)=0,          rhs_y(nx+2,ny+1,nz+2)=0,          rhs_z(nx+2,ny+2,nz+1)=0
    !real(8) :: rhs_x_previous(nx+1,ny+2,nz+2)=0, rhs_y_previous(nx+2,ny+1,nz+2)=0, rhs_z_previous(nx+2,ny+2,nz+1)=0
    !real(8) :: dpdx(nx+1,ny+2,nz+2)=0,           dpdy(nx+2,ny+1,nz+2)=0,           dpdz(nx+2,ny+2,nz+1)=0
    real(8) :: f_term_x=0,                       f_term_y=0,                       f_term_z=0
    !real(8) :: rhs_x_previous0(nx+1,ny+2,nz+2)=0,rhs_y_previous0(nx+2,ny+1,nz+2)=0,rhs_z_previous0(nx+2,ny+2,nz+1)=0
    !real(8) :: f_term_x(nx+1,ny+2,nz+2)=0,       f_term_y(nx+2,ny+1,nz+2)=0,       f_term_z(nx+2,ny+2,nz+1)=0
    real(8), dimension (:,:,:), allocatable :: rhs_x, rhs_y, rhs_z, rhs_x_previous, rhs_y_previous, rhs_z_previous, dpdx, dpdy, dpdz
    real(8), dimension (:,:,:), allocatable :: dp_lu
    real(8), dimension (:    ), allocatable :: RHS_poisson0, dp_vec
    real(8), dimension (:,:,:), allocatable :: conv_x, conv_y, conv_z
    real(8), dimension (:,:,:), allocatable :: diff_x, diff_y, diff_z
    real(8), dimension (:,:,:), allocatable :: RHS_poisson_internal, div

    !!!!!!!!!!!!!!! for CN2-ADI scheme !!!!!!!!!!!!!!!
    real(8), dimension (:), allocatable :: A_low, A_d, A_up, A_up2, B_low, B_d, B_up, B_up2, C_low, C_d, C_up, C_up2
    real(8), dimension (:,:), allocatable :: A_mat, B_mat, C_mat, D_mat, E_mat, F_mat
    integer, dimension (:), allocatable :: A_ipiv, B_ipiv, C_ipiv, D_ipiv, E_ipiv, F_ipiv

    !!!!!!!!!!!!!!! for re-simulation !!!!!!!!!!!!!!!
    real(8) :: u_sub(nx0+1,ny0+2,nz0+2)=0,       v_sub(nx0+2,ny0+1,nz0+2)=0,       w_sub(nx0+2,ny0+2,nz0+1)=0,       p_sub(nx0+2,ny0+2,nz0+2)=0
    real(8) :: u_star_sub(nx0+1,ny0+2,nz0+2)=0,  v_star_sub(nx0+2,ny0+1,nz0+2)=0,  w_star_sub(nx0+2,ny0+2,nz0+1)=0,  dp_sub(nx0+2,ny0+2,nz0+2)=0
    real(8) :: RHS_poisson_sub(nx0+2,ny0+2,nz0+2)=0

    !!!!!!!!!!!!!!! boundary conditions !!!!!!!!!!!!!!!
    real(8) :: bx_u_1(ny+2,nz+2)=0, bx_u_nx(ny+2,nz+2)=0, by_u_1(nx+1,nz+2)=0, by_u_ny(nx+1,nz+2)=0, bz_u_1(nx+1,ny+2)=0, bz_u_nz(nx+1,ny+2)=0
    real(8) :: bx_v_1(ny+1,nz+2)=0, bx_v_nx(ny+1,nz+2)=0, by_v_1(nx+2,nz+2)=0, by_v_ny(nx+2,nz+2)=0, bz_v_1(nx+2,ny+1)=0, bz_v_nz(nx+2,ny+1)=0
    real(8) :: bx_w_1(ny+2,nz+1)=0, bx_w_nx(ny+2,nz+1)=0, by_w_1(nx+2,nz+1)=0, by_w_ny(nx+2,nz+1)=0, bz_w_1(nx+2,ny+2)=0, bz_w_nz(nx+2,ny+2)=0
    real(8) :: bx_p_1(nyp,nzp)=0,   bx_p_nx(nyp,nzp)=0,   by_p_1(nxp,nzp)=0,   by_p_ny(nxp,nzp)=0,   bz_p_1(nxp,nyp)=0,   bz_p_nz(nxp,nyp)=0

    real(8) :: bx_u_1s(ny0+2,nz0+2)=0, bx_u_nxs(ny0+2,nz0+2)=0, by_u_1s(nx0+1,nz0+2)=0, by_u_nys(nx0+1,nz0+2)=0, bz_u_1s(nx0+1,ny0+2)=0, bz_u_nzs(nx0+1,ny0+2)=0
    real(8) :: bx_v_1s(ny0+1,nz0+2)=0, bx_v_nxs(ny0+1,nz0+2)=0, by_v_1s(nx0+2,nz0+2)=0, by_v_nys(nx0+2,nz0+2)=0, bz_v_1s(nx0+2,ny0+1)=0, bz_v_nzs(nx0+2,ny0+1)=0
    real(8) :: bx_w_1s(ny0+2,nz0+1)=0, bx_w_nxs(ny0+2,nz0+1)=0, by_w_1s(nx0+2,nz0+1)=0, by_w_nys(nx0+2,nz0+1)=0, bz_w_1s(nx0+2,ny0+2)=0, bz_w_nzs(nx0+2,ny0+2)=0

    real(8) :: tGet
    integer :: time_length, t_step, plot_step=20, slice=nz/2+1
    real(8), dimension (:), allocatable :: time_array
    integer, dimension (:), allocatable :: x_range0, y_range0, z_range0
    integer :: x_range1(nx)=[(i, i=2, nx+1)], y_range1(ny)=[(i, i=2, ny+1)], z_range1(nz)=[(i, i=2, nz+1)]
    type(csr), allocatable :: LHS_poisson
    !type(coo), allocatable :: LHS_poisson_coo
    !real(8), dimension (:,:), allocatable :: LHS_poisson_den

    !!!!!!!!!!!!!!! temperary variables !!!!!!!!!!!!!!!
    !integer(8) :: sizeof_record, sizeof_record_sub, sizeof_record_slice, tempi1, tempi2
    !real(8) :: temp01, temp02, temp03, temp04, temp05, temp06
    !real(8), dimension (:), allocatable :: temp11, temp12, temp13, temp14, temp15, temp16
    real(8), dimension (:,:), allocatable :: temp21!, temp22, temp23, temp24, temp25, temp26
    !real(8), dimension (:,:,:), allocatable :: temp31, temp32, temp33, temp34, temp35, temp36
    !real(4) :: ini_vel(3,nx_file,nx_file,nx_file), ini_pr(nx_file,nx_file,nx_file)
    !real(4), dimension (:), allocatable :: temp11s
    !real(4) :: temp01s

    !!!!!!!!!!!!!!! INTEL mkl_pardiso !!!!!!!!!!!!!!!
    type(MKL_PARDISO_HANDLE) pt(64)
    integer :: maxfct=1, mnum=1, mtype=11, phase=13, n=nx*ny*nz, idum(nx*ny*nz), nrhs=1, iparm(64)=0, msglvl=0, error=0

    !!!!!!!!!!!!!!! INTEL mkl_dft !!!!!!!!!!!!!!!
    integer :: cstrides(4)=0, rstrides(4)=0
    type(DFTI_DESCRIPTOR), POINTER :: hand_f, hand_b
    !real(8) :: dft_out_r(nx,ny,nz), poisson_eigv(nx,ny,nz)=0
    !complex(8) :: dft_out_c(nx,ny,nz)
    real(8) :: dft_out_r1(nx,ny,nz), poisson_eigv(INT(nx/2.0)+1,ny,nz)=0
    complex(8) :: dft_out_c1(INT(nx/2.0)+1,ny,nz)
    integer :: status

    !!!!!!!!!!!!!!! HDF5 !!!!!!!!!!!!!!!
    character(len=1000) :: string_var
    integer(8) :: h5f_whole, h5f_sub, h5f_slice, h5g_sub, h5g_slice, prp_id
    logical :: first_write=.true.

    !!!!!!!!!!!!!!! system_clock !!!!!!!!!!!!!!!
    REAL(8) :: system_clock_rate
    INTEGER :: c01,c02,c1,c2,cr,cm

    !!!!!!!!!!!!!!! simulation parameters !!!!!!!!!!!!!!!
    real(8) :: dt0=0.004d0
    character(*), parameter :: timescheme="AB2"
    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall); pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    integer, parameter :: bc_x=1, bc_y=bc_x, bc_z=bc_x, pbc_x=1, pbc_y=pbc_x, pbc_z=pbc_x
    logical, parameter :: save_output=.true., LU_poisson=(nxp*nyp*nzp<=34**3), FFT_poisson=(pbc_x==1 .and. pbc_y==1 .and. pbc_z==1)

    call OMP_set_dynamic(.true.)
    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)
    ! open HDF5 module
    call h5open_f(status)

    call get_command_argument (1, string_var, status=status)
    if (status==0) then
        read(string_var,*) dt0
    end if
    time_length=nint((t_end-t_start)/(dt0))
    allocate( time_array(0:time_length) )
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

    allocate( rhs_x(nx+1,ny+2,nz+2),          rhs_y(nx+2,ny+1,nz+2),          rhs_z(nx+2,ny+2,nz+1) )
    allocate( rhs_x_previous(nx+1,ny+2,nz+2), rhs_y_previous(nx+2,ny+1,nz+2), rhs_z_previous(nx+2,ny+2,nz+1) )
    allocate( dpdx(nx+1,ny+2,nz+2),           dpdy(nx+2,ny+1,nz+2),           dpdz(nx+2,ny+2,nz+1))
    allocate( div(nx,ny,nz), RHS_poisson_internal(nx,ny,nz) )

    if (timescheme=="AB2-CN") then
        if (allocated(A_mat)) deallocate(A_mat, B_mat, C_mat, D_mat, E_mat, F_mat)
        if (allocated(A_ipiv)) deallocate(A_ipiv, B_ipiv, C_ipiv, D_ipiv, E_ipiv, F_ipiv)
        allocate(A_mat(nx+1,nx+1), B_mat(nx+2,nx+2), C_mat(ny+1,ny+1), D_mat(ny+2,ny+2), E_mat(nz+1,nz+1), F_mat(nz+2,nz+2))
        allocate(A_ipiv(nx+1), B_ipiv(nx+2), C_ipiv(ny+1), D_ipiv(ny+2), E_ipiv(nz+1), F_ipiv(nz+2))

        call mat_CN(nx, bc_x, dx2, dx, A_mat, B_mat, dt0, nu)
        call mat_CN(ny, bc_y, dy2, dy, C_mat, D_mat, dt0, nu)
        call mat_CN(nz, bc_z, dz2, dz, E_mat, F_mat, dt0, nu)

        call getrf( A_mat, A_ipiv )
        call getrf( B_mat, B_ipiv )
        call getrf( C_mat, C_ipiv )
        call getrf( D_mat, D_ipiv )
        call getrf( E_mat, E_ipiv )
        call getrf( F_mat, F_ipiv )
    end if

    if (LU_poisson) then
        allocate(RHS_poisson0(nx*ny*nz), dp_vec(nx*ny*nz), dp_lu(nxp,nyp,nzp))

        LHS_poisson=Poisson_LHS_staggered(nx, ny, nz, dx2, dy2, dz2, pbc_x, pbc_y, pbc_z, dx, dy, dz)
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
        poisson_eigv(1,1,1)=1.0d0
        poisson_eigv=-4.0d0*poisson_eigv
    end if

    do t_step=0,time_length
        tGet=time_array(t_step)
        print *,''
        print '(" t_step ", I6, "        tGet ", F7.4)', t_step, tGet
        CALL SYSTEM_CLOCK(c01)

        if (t_step==0) then
            if (init) then
                !if (allocated(temp11s)) deallocate(temp11s)
                !allocate(temp11s(size(u)+size(v)+size(w)+size(p)))
                !write (string_var,'("D:\Documents\source\repos\NS_3D_Staggered_dp\NS_3D_Staggered_dp\256^3_init_from_JHTDB.dat")')
                !open(20, file=string_var, form='unformatted', status='old', action='read', &
                !    access='direct', recl=size(u)+size(v)+size(w)+size(p)) !variables are single precision
                !read(20, rec=1) temp11s
                !close(20)
                !
                !tempi1=1;        tempi2=size(u);        u=reshape(dble(temp11s(tempi1:tempi2)), [nx+1,ny+2,nz+2]);
                !tempi1=tempi2+1; tempi2=tempi2+size(v); v=reshape(dble(temp11s(tempi1:tempi2)), [nx+2,ny+1,nz+2]);
                !tempi1=tempi2+1; tempi2=tempi2+size(w); w=reshape(dble(temp11s(tempi1:tempi2)), [nx+2,ny+2,nz+1]);
                !tempi1=tempi2+1; tempi2=tempi2+size(p); p=reshape(dble(temp11s(tempi1:tempi2)), [nx+2,ny+2,nz+2]);
                !
                !deallocate(temp11s)

                !sizeof_record=size(u)+size(v)+size(w)+size(p)+size(u)+size(v)+size(w)
                !
                !if (allocated(temp11)) deallocate(temp11)
                !allocate(temp11(sizeof_record))
                !write (string_var,'("D:\Documents\source\repos\NS_3D_Staggered_dp\NS_3D_Staggered_dp\init_file\HIT_256^3_decay_4.E-3_AB2_dp_init.dat")')
                !open(20, file=string_var, form='unformatted', status='old', action='read', &
                !    access='direct', recl=sizeof_record*2) !variables are double precision
                !read(20, rec=1) temp11
                !close(20)
                !
                !tempi1=1;        tempi2=size(u);        u_star=reshape(temp11(tempi1:tempi2), [nx+1,ny+2,nz+2]);
                !tempi1=tempi2+1; tempi2=tempi2+size(v); v_star=reshape(temp11(tempi1:tempi2), [nx+2,ny+1,nz+2]);
                !tempi1=tempi2+1; tempi2=tempi2+size(w); w_star=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+1]);
                !tempi1=tempi2+1; tempi2=tempi2+size(p); p=reshape(temp11(tempi1:tempi2), [nx+2,ny+2,nz+2]);

                call h5fopen_f("init_file/HIT_256^3_decay_4.E-3_AB2_dp_init.h5", H5F_ACC_RDONLY_F, h5f_whole, status)

                call h5ltread_dataset_double_f(h5f_whole, 'u', u, [nx+1_8,ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'v', v, [nx+2_8,ny+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'w', w, [nx+2_8,ny+2_8,nz+1_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'p', p, [nx+2_8,ny+2_8,nz+2_8], status)

                call h5fclose_f(h5f_whole, status)

                !call TGV(xu, yu, zu, 0.0d0, nu, u)
                !call TGV(xv, yv, zv, 0.0d0, nu, v=v)
                !call TGV(xp, yp, zp, 0.0d0, nu, p=p)
            else
                write (string_var,'(A, "_result/HIT_", I0, "^3_decay_", ES5.0E1, "_", A ,"_dp_t_" , F0.4, ".h5")') trim(timescheme), nx, dt0, trim(timescheme), tGet
                call h5fopen_f(string_var, H5F_ACC_RDONLY_F, h5f_whole, status)

                call h5ltread_dataset_double_f(h5f_whole, 'u', u, [nx+1_8,ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'v', v, [nx+2_8,ny+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'w', w, [nx+2_8,ny+2_8,nz+1_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'p', p, [nx+2_8,ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'rhs_x_previous', rhs_x_previous, [nx+1_8,ny+2_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'rhs_y_previous', rhs_y_previous, [nx+2_8,ny+1_8,nz+2_8], status)
                call h5ltread_dataset_double_f(h5f_whole, 'rhs_z_previous', rhs_z_previous, [nx+2_8,ny+2_8,nz+1_8], status)

                call h5fclose_f(h5f_whole, status)
            end if
        else
            !!! Convection
            CALL SYSTEM_CLOCK(c1)
            call cal_conv(u, v, w, bc_x, bc_y, bc_z, dx, dy, dz, conv_x, conv_y, conv_z)
            CALL SYSTEM_CLOCK(c2)
            print '("    convective term completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

            !!! Diffusion
            CALL SYSTEM_CLOCK(c1)
            call cal_diff(u, v, w, bc_x, bc_y, bc_z, dx2, dy2, dz2, diff_x, diff_y, diff_z)
            CALL SYSTEM_CLOCK(c2)
            print '("    diffusion term completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

            !!! Time-advancement
            f_term_x=0; f_term_y=0; f_term_z=0;

            CALL SYSTEM_CLOCK(c1)
            dpdx=diff(p,1,1)/dx;
            dpdy=diff(p,1,2)/dy;
            dpdz=diff(p,1,3)/dz;
            CALL SYSTEM_CLOCK(c2)
            print '("    pressure gradient completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

            CALL SYSTEM_CLOCK(c1)
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
                !if (t_step==1 .and. init) then
                if ( (t_step==1 .or. t_step==nint( (dt0-t_start)/dt0) ) .and. init) then
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

                rhs_x_previous(x_range0,y_range1,z_range1)=-conv_x;
                rhs_y_previous(x_range1,y_range0,z_range1)=-conv_y;
                rhs_z_previous(x_range1,y_range1,z_range0)=-conv_z;

                call vel_bc_staggered_CN2(u_star, v_star, w_star, bx_u_1, bx_u_nx, bx_v_1, bx_v_nx, bx_w_1, bx_w_nx, bc_x, 1)
                !!$omp parallel do
                do k=1,nz+2
                    call getrs( A_mat, A_ipiv, u_star(:,:,k) )
                    call getrs( B_mat, B_ipiv, v_star(:,:,k) )
                    if (k<=nz+1) call getrs( B_mat, B_ipiv, w_star(:,:,k) )
                end do
                !!$omp end parallel do
                bx_u_1s = u_star(idx_xu(1),     idx_yu, idx_zu)
                bx_u_nxs= u_star(idx_xu(nx0+1), idx_yu, idx_zu)
                bx_v_1s =(v_star(idx_xv(1),     idx_yv, idx_zv) + v_star(idx_xv(2),     idx_yv, idx_zv))*0.5d0
                bx_v_nxs=(v_star(idx_xv(nx0+1), idx_yv, idx_zv) + v_star(idx_xv(nx0+2), idx_yv, idx_zv))*0.5d0
                bx_w_1s =(w_star(idx_xw(1),     idx_yw, idx_zw) + w_star(idx_xw(2),     idx_yw, idx_zw))*0.5d0
                bx_w_nxs=(w_star(idx_xw(nx0+1), idx_yw, idx_zw) + w_star(idx_xw(nx0+2), idx_yw, idx_zw))*0.5d0
                !temp31=u_star(idx_xu,idx_yu,idx_zu); temp32=v_star(idx_xv,idx_yv,idx_zv); temp33=w_star(idx_xw,idx_yw,idx_zw);

                call vel_bc_staggered_CN2(u_star, v_star, w_star, by_u_1, by_u_ny, by_v_1, by_v_ny, by_w_1, by_w_ny, bc_y, 2)
                !!$omp parallel do
                do k=1,nz+2
                    temp21=transpose(u_star(:,:,k))
                    call getrs( D_mat, D_ipiv, temp21 )
                    u_star(:,:,k)=transpose(temp21)

                    temp21=transpose(v_star(:,:,k))
                    call getrs( C_mat, C_ipiv, temp21 )
                    v_star(:,:,k)=transpose(temp21)

                    if (k<=nz+1) then
                        temp21=transpose(w_star(:,:,k))
                        call getrs( D_mat, D_ipiv, temp21 )
                        w_star(:,:,k)=transpose(temp21)
                    end if
                end do
                !!$omp end parallel do
                by_u_1s =(u_star(idx_xu, idx_yu(1),     idx_zu) + u_star(idx_xu, idx_yu(2),     idx_zu))*0.5d0
                by_u_nys=(u_star(idx_xu, idx_yu(ny0+1), idx_zu) + u_star(idx_xu, idx_yu(ny0+2), idx_zu))*0.5d0
                by_v_1s = v_star(idx_xv, idx_yv(1),     idx_zv)
                by_v_nys= v_star(idx_xv, idx_yv(ny0+1), idx_zv)
                by_w_1s =(w_star(idx_xw, idx_yw(1),     idx_zw) + w_star(idx_xw, idx_yw(2),     idx_zw))*0.5d0
                by_w_nys=(w_star(idx_xw, idx_yw(ny0+1), idx_zw) + w_star(idx_xw, idx_yw(ny0+2), idx_zw))*0.5d0

                call vel_bc_staggered_CN2(u_star, v_star, w_star, bz_u_1, bz_u_nz, bz_v_1, bz_v_nz, bz_w_1, bz_w_nz, bc_z, 3)
                !!$omp parallel do
                do j=1,ny+2
                    temp21=transpose(u_star(:,j,:))
                    call getrs( F_mat, F_ipiv, temp21 )
                    u_star(:,j,:)=transpose(temp21)

                    if (j<=ny+1) then
                        temp21=transpose(v_star(:,j,:))
                        call getrs( F_mat, F_ipiv, temp21 )
                        v_star(:,j,:)=transpose(temp21)
                    end if

                    temp21=transpose(w_star(:,j,:))
                    call getrs( E_mat, E_ipiv, temp21 )
                    w_star(:,j,:)=transpose(temp21)
                end do
                !!$omp end parallel do
                bz_u_1s =(u_star(idx_xu, idx_yu, idx_zu(1))     + u_star(idx_xu, idx_yu, idx_zu(2))    )*0.5d0
                bz_u_nzs=(u_star(idx_xu, idx_yu, idx_zu(nz0+1)) + u_star(idx_xu, idx_yu, idx_zu(nz0+2)))*0.5d0
                bz_v_1s =(v_star(idx_xv, idx_yv, idx_zv(1))     + v_star(idx_xv, idx_yv, idx_zv(2))    )*0.5d0
                bz_v_nzs=(v_star(idx_xv, idx_yv, idx_zv(nz0+1)) + v_star(idx_xv, idx_yv, idx_zv(nz0+2)))*0.5d0
                bz_w_1s = w_star(idx_xw, idx_yw, idx_zw(1))
                bz_w_nzs= w_star(idx_xw, idx_yw, idx_zw(nz0+1))
            end if
            CALL SYSTEM_CLOCK(c2)
            print '("    time-advancement completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

            !RHS_poisson(2:ubound(RHS_poisson,1)-1,2:ubound(RHS_poisson,2)-1,2:ubound(RHS_poisson,3)-1) = &
            !    (diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx + &
            !    diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy + &
            !    diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz) / (dt0);
            RHS_poisson_internal = diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx
            RHS_poisson_internal = RHS_poisson_internal + diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy
            RHS_poisson_internal = RHS_poisson_internal + diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz
            RHS_poisson_internal = RHS_poisson_internal/dt0

            call pr_bc_staggered_modify_rhs(RHS_poisson_internal, bx_p_1(2:nyp-1,2:nzp-1), bx_p_nx(2:nyp-1,2:nzp-1), &
                by_p_1(2:nxp-1,2:nzp-1), by_p_ny(2:nxp-1,2:nzp-1), bz_p_1(2:nxp-1,2:nyp-1), bz_p_nz(2:nxp-1,2:nyp-1), &
                pbc_x, pbc_y, pbc_z, dx, dy, dz, dx2, dy2, dz2)

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
                dft_out_c1(1,1,1)=0
                !print *, "   Backward FFT..."
                status = DftiComputeBackward(hand_b, dft_out_c1(:,1,1), dft_out_r1(:,1,1))
                !status = DftiComputeBackward(hand_b, dft_out_c(:,1,1))
                !dft_out_r=dble(dft_out_c)
                dp(2:nxp-1,2:nyp-1,2:nzp-1)=dft_out_r1
                call pr_bc_staggered(dp, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz, pbc_x, pbc_y, pbc_z, dx, dy, dz)

                CALL SYSTEM_CLOCK(c2)
                print '("    Solve Poisson (FFT-based FD) completed: ", F8.4, " second")', (c2-c1)/system_clock_rate
                !print *, "**************************************"
            end if

            !print *, maxval(abs(dp-dp(1,1,1)-dp_lu+dp_lu(1,1,1)))
            !temp31=dp-dp_lu

            !OPEN(10, file="poisson_eq.dat", form="unformatted")
            !WRITE(10) [RHS_poisson]
            !WRITE(10) [dp]
            !WRITE(10) [dp_lu]
            !CLOSE(10)
            !dp=LHS_poisson\RHS_poisson(:);
            !dp=reshape(dp,nxp,nyp,nzp);

            !if (LU_poisson) dp=dp_lu
            dp=dp-dp(1,1,1)
            CALL SYSTEM_CLOCK(c1)
            u=-dt0*diff(dp,1,1)/dx+u_star
            v=-dt0*diff(dp,1,2)/dy+v_star
            w=-dt0*diff(dp,1,3)/dz+w_star
            p=p+dp
            CALL SYSTEM_CLOCK(c2)
            print '("    correction completed: ", F8.4, " second")', (c2-c1)/system_clock_rate

        end if
        print *, u(nx/2,ny/2,nz/2)
        CALL SYSTEM_CLOCK(c02)

        !!!!!!!!!!!!!!!!!!!!!!!!!! post-processing !!!!!!!!!!!!!!!!!!!!!!!!!!
        div=diff_old(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,1)/dx + &
            diff_old(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),1,2)/dy + &
            diff_old(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),1,3)/dz;
        print '("Complete: ", F8.4, " second. MAX Div: ", E13.6)', (c02-c01)/system_clock_rate, maxval(abs(div))
        print *, "**************************************"

        if (tGet>=0 .and. save_output) then

            if (first_write) then
                !call h5pcreate_f(H5P_FILE_ACCESS_F, prp_id, status)
                !call h5pset_fapl_core_f(prp_id, 64, .true., status)
                !call h5pset_cache_f(prp_id, 0, 521, 16000000, 1.0, status)

                write (string_var,'(A, "_result/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_dp_x0_", I0, "_nx0_", I0, "_sub.h5")') trim(timescheme), nx, dt0, trim(timescheme), x0, nx0
                call h5fcreate_f(string_var, H5F_ACC_EXCL_F, h5f_sub, status)
                call h5ltset_attribute_int_f(h5f_sub, "/", "nx", [nx], 1, status)
                call h5ltset_attribute_double_f(h5f_sub, "/", "dt", [dt0], 1, status)
                call h5ltset_attribute_double_f(h5f_sub, "/", "nu", [nu], 1, status)
                call h5ltset_attribute_string_f(h5f_sub, "/", "time_scheme", trim(timescheme), status)
                call h5ltset_attribute_int_f(h5f_sub, "/", "time_step", [t_step], 1, status)
                call h5ltset_attribute_double_f(h5f_sub, "/", "time", [tGet], 1, status)
                call h5ltset_attribute_int_f(h5f_sub, "/", "x0", [x0], 1, status)
                call h5ltset_attribute_int_f(h5f_sub, "/", "nx0", [nx0], 1, status)

                write (string_var,'(A, "_result/HIT_", I0,"^3_decay_", ES5.0E1, "_", A , "_slice_", I0, ".h5")') trim(timescheme), nx, dt0, trim(timescheme), slice
                call h5fcreate_f(string_var, H5F_ACC_EXCL_F, h5f_slice, status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "nx", [nx], 1, status)
                call h5ltset_attribute_double_f(h5f_slice, "/", "dt", [dt0], 1, status)
                call h5ltset_attribute_double_f(h5f_slice, "/", "nu", [nu], 1, status)
                call h5ltset_attribute_string_f(h5f_slice, "/", "time_scheme", trim(timescheme), status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "time_step", [t_step], 1, status)
                call h5ltset_attribute_double_f(h5f_slice, "/", "time", [tGet], 1, status)
                call h5ltset_attribute_int_f(h5f_slice, "/", "XY-slice", [slice], 1, status)
                
                first_write=.false.
            end if

            if (mod(abs(tGet), 1.0d0) < dt0/sub_tstep/100.0d0) then
                if (init .or. (.not. init .and. t_step/=0)) then
                    write (string_var,'(A, "_result/HIT_", I0, "^3_decay_", ES5.0E1, "_", A ,"_dp_t_" , F0.4, ".h5")') trim(timescheme), nx, dt0, trim(timescheme), tGet
                    call h5fcreate_f(string_var, H5F_ACC_EXCL_F, h5f_whole, status)

                    call h5ltset_attribute_int_f(h5f_whole, "/", "nx", [nx], 1, status)
                    call h5ltset_attribute_double_f(h5f_whole, "/", "dt", [dt0], 1, status)
                    call h5ltset_attribute_double_f(h5f_whole, "/", "nu", [nu], 1, status)
                    call h5ltset_attribute_string_f(h5f_whole, "/", "time_scheme", trim(timescheme), status)
                    call h5ltset_attribute_int_f(h5f_whole, "/", "time_step", [t_step], 1, status)
                    call h5ltset_attribute_double_f(h5f_whole, "/", "time", [tGet], 1, status)

                    call h5ltmake_dataset_double_f(h5f_whole, 'u', 3, [nx+1_8,ny+2_8,nz+2_8], u, status)
                    call h5ltmake_dataset_double_f(h5f_whole, 'v', 3, [nx+2_8,ny+1_8,nz+2_8], v, status)
                    call h5ltmake_dataset_double_f(h5f_whole, 'w', 3, [nx+2_8,ny+2_8,nz+1_8], w, status)
                    call h5ltmake_dataset_double_f(h5f_whole, 'p', 3, [nx+2_8,ny+2_8,nz+2_8], p, status)
                    call h5ltmake_dataset_double_f(h5f_whole, 'rhs_x_previous', 3, [nx+1_8,ny+2_8,nz+2_8], rhs_x_previous, status)
                    call h5ltmake_dataset_double_f(h5f_whole, 'rhs_y_previous', 3, [nx+2_8,ny+1_8,nz+2_8], rhs_y_previous, status)
                    call h5ltmake_dataset_double_f(h5f_whole, 'rhs_z_previous', 3, [nx+2_8,ny+2_8,nz+1_8], rhs_z_previous, status)

                    call h5fclose_f(h5f_whole, status)

                    call h5fclose_f(h5f_sub, status)
                    write (string_var,'(A, "_result/HIT_", I0, "^3_decay_", ES5.0E1, "_", A , "_dp_x0_", I0, "_nx0_", I0, "_sub.h5")') trim(timescheme), nx, dt0, trim(timescheme), x0, nx0
                    call h5fopen_f(string_var, H5F_ACC_RDWR_F, h5f_sub, status)

                    call h5fclose_f(h5f_slice, status)
                    write (string_var,'(A, "_result/HIT_", I0,"^3_decay_", ES5.0E1, "_", A , "_slice_", I0, ".h5")') trim(timescheme), nx, dt0, trim(timescheme), slice
                    call h5fopen_f(string_var, H5F_ACC_RDWR_F, h5f_slice, status)
                end if
            end if

            if (tGet<=20) then
                u_sub=u(idx_xu,idx_yu,idx_zu); u_star_sub=u_star(idx_xu,idx_yu,idx_zu);
                v_sub=v(idx_xv,idx_yv,idx_zv); v_star_sub=v_star(idx_xv,idx_yv,idx_zv);
                w_sub=w(idx_xw,idx_yw,idx_zw); w_star_sub=w_star(idx_xw,idx_yw,idx_zw);
                p_sub=p(idx_xp,idx_yp,idx_zp); dp_sub=dp(idx_xp,idx_yp,idx_zp)
                RHS_poisson_sub=RHS_poisson_internal(idx_xp-1,idx_yp-1,idx_zp-1) !RHS_poisson(idx_xp,idx_yp,idx_zp)

                write (string_var,'("t_", F0.4)') tGet
                call h5gcreate_f(h5f_sub, string_var, h5g_sub, status)

                call h5ltmake_dataset_double_f(h5g_sub, 'u_sub', 3, [nx0+1_8,ny0+2_8,nz0+2_8], u_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'v_sub', 3, [nx0+2_8,ny0+1_8,nz0+2_8], v_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'w_sub', 3, [nx0+2_8,ny0+2_8,nz0+1_8], w_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'p_sub', 3, [nx0+2_8,ny0+2_8,nz0+2_8], p_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'u_star_sub', 3, [nx0+1_8,ny0+2_8,nz0+2_8], u_star_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'v_star_sub', 3, [nx0+2_8,ny0+1_8,nz0+2_8], v_star_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'w_star_sub', 3, [nx0+2_8,ny0+2_8,nz0+1_8], w_star_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'dp_sub', 3, [nx0+2_8,ny0+2_8,nz0+2_8], dp_sub, status)
                call h5ltmake_dataset_double_f(h5g_sub, 'RHS_poisson_sub', 3, [nx0+2_8,ny0+2_8,nz0+2_8], RHS_poisson_sub, status)

                if (mod(tGet,1.0d0)==0.0d0) then
                    call h5ltmake_dataset_double_f(h5g_sub, 'rhs_x_previous', 3, [nx0+1_8,ny0+2_8,nz0+2_8], rhs_x_previous(idx_xu,idx_yu,idx_zu), status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'rhs_y_previous', 3, [nx0+2_8,ny0+1_8,nz0+2_8], rhs_y_previous(idx_xv,idx_yv,idx_zv), status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'rhs_z_previous', 3, [nx0+2_8,ny0+2_8,nz0+1_8], rhs_z_previous(idx_xw,idx_yw,idx_zw), status)
                end if
                if (timescheme=="AB2-CN") then
                    call h5ltmake_dataset_double_f(h5g_sub, 'bx_u_1s', 2, [ny0+2_8,nz0+2_8], bx_u_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bx_u_nxs', 2, [ny0+2_8,nz0+2_8], bx_u_nxs, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bx_v_1s', 2, [ny0+1_8,nz0+2_8], bx_v_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bx_v_nxs', 2, [ny0+1_8,nz0+2_8], bx_v_nxs, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bx_w_1s', 2, [ny0+2_8,nz0+1_8], bx_w_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bx_w_nxs', 2, [ny0+2_8,nz0+1_8], bx_w_nxs, status)

                    call h5ltmake_dataset_double_f(h5g_sub, 'by_u_1s', 2, [nx0+1_8,nz0+2_8], by_u_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'by_u_nys', 2, [nx0+1_8,nz0+2_8], by_u_nys, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'by_v_1s', 2, [nx0+2_8,nz0+2_8], by_v_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'by_v_nys', 2, [nx0+2_8,nz0+2_8], by_v_nys, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'by_w_1s', 2, [nx0+2_8,nz0+1_8], by_w_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'by_w_nys', 2, [nx0+2_8,nz0+1_8], by_w_nys, status)

                    call h5ltmake_dataset_double_f(h5g_sub, 'bz_u_1s', 2, [nx0+1_8,ny0+2_8], bz_u_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bz_u_nzs', 2, [nx0+1_8,ny0+2_8], bz_u_nzs, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bz_v_1s', 2, [nx0+2_8,ny0+1_8], bz_v_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bz_v_nzs', 2, [nx0+2_8,ny0+1_8], bz_v_nzs, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bz_w_1s', 2, [nx0+2_8,ny0+2_8], bz_w_1s, status)
                    call h5ltmake_dataset_double_f(h5g_sub, 'bz_w_nzs', 2, [nx0+2_8,ny0+2_8], bz_w_nzs, status)
                end if
                call h5gclose_f( h5g_sub, status)
            end if

            if (mod(t_step,plot_step)==0 .and. tGet<=20) then
                write (string_var,'("t_", F0.4)') tGet
                call h5gcreate_f(h5f_slice, string_var, h5g_slice, status)

                call h5ltmake_dataset_double_f(h5g_slice, 'u_slice', 2, [nx+1_8,ny+2_8], u(:,:,slice), status)
                call h5ltmake_dataset_double_f(h5g_slice, 'v_slice', 2, [nx+2_8,ny+1_8], v(:,:,slice), status)
                call h5ltmake_dataset_double_f(h5g_slice, 'w_slice', 2, [nx+2_8,ny+2_8], w(:,:,slice), status)
                call h5ltmake_dataset_double_f(h5g_slice, 'p_slice', 2, [nx+2_8,ny+2_8], p(:,:,slice), status)

                call h5gclose_f( h5g_slice, status)
            end if

        end if


    end do

    if (FFT_poisson) then
        status = DftiFreeDescriptor(hand_f)
        status = DftiFreeDescriptor(hand_b)
    end if

    call h5fclose_f(h5f_sub, status)
    call h5fclose_f(h5f_slice, status)
    call h5close_f(status)

    end subroutine big_periodic_3D