    !  NS_3D_Staggered_dp.f90
    !
    !  FUNCTIONS:
    !  NS_3D_Staggered_dp - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: NS_3D_Staggered_dp
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************
    include 'mkl_dfti.f90'

    program NS_3D_Staggered_dp

    use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R
    use FD_functions
    use NS_functions
    use coo_mod
    use csr_mod
    !use sparse_pack
    use other_utility
    use ogpf

    implicit none
    include 'omp_lib.h'
    INCLUDE 'mkl_pardiso.fi'
    include 'fftw/fftw3.f'

    ! Variables
    real(8), parameter :: pi = 3.1415926535897932_8
    integer :: i=0,j=0,k=0,ll=0,mm=0

    logical, parameter :: init=.true.
    character *4 :: timescheme="AB2"
    !integer, parameter :: timescheme=2 ! 1 Euler; 2 AB2
    real(8), parameter :: Re=1.0d0, nu=0.005d0, t_start=-1.0d0, t_end=10.0d0, dt0=0.005d0
    integer, parameter :: nx=128, ny=nx, nz=nx, nxp=nx+2, nyp=ny+2, nzp=nz+2, sub_tstep=1
    logical, parameter :: LU_poisson=(nxp*nyp*nzp<=34**3)
    real(8), parameter :: lx=2.0d0*pi, ly=2.0d0*pi, lz=2.0d0*pi, dx=lx/nx, dy=ly/ny, dz=lz/nz, dx2=dx*dx, dy2=dy*dy, dz2=dz*dz
    real(8), parameter :: xu(nx+1)=[(i*dx, i=0, nx)],        yu(ny+2)=[((i+0.5)*dy, i=-1, ny)], zu(nz+2)=[((i+0.5)*dz, i=-1, nz)]
    real(8), parameter :: xv(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yv(ny+1)=[(i*dy, i=0, ny)],          zv(nz+2)=[((i+0.5d0)*dz, i=-1, nz)]
    real(8), parameter :: xw(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yw(ny+2)=[((i+0.5d0)*dy, i=-1, ny)], zw(nz+1)=[(i*dz, i=0, nz)]
    real(8), parameter :: xp(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yp(ny+2)=[((i+0.5d0)*dy, i=-1, ny)], zp(nz+2)=[((i+0.5d0)*dz, i=-1, nz)]

    !integer, parameter :: nx0=8, ny0=8, nz0=8, x0=9, y0=9, z0=9
    !integer, dimension (:), allocatable :: idx_xu, idx_yu, idx_zu, idx_xv, idx_yv, idx_zv, idx_xw, idx_yw, idx_zw, idx_xp, idx_yp, idx_zp

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall)
    ! pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    integer, parameter :: bc_x=1, bc_y=1, bc_z=1, pbc_x=1, pbc_y=1, pbc_z=1

    integer, parameter :: time_length=nint((t_end-t_start)/(dt0))
    real(8), dimension (0:time_length), parameter :: time_array=[(t_start+dt0*i, i=0, time_length)]
    real(8) :: u(nx+1,ny+2,nz+2)=0,              v(nx+2,ny+1,nz+2)=0,              w(nx+2,ny+2,nz+1)=0,       p(nxp,nyp,nzp)=0
    real(8) :: u_star(nx+1,ny+2,nz+2)=0,         v_star(nx+2,ny+1,nz+2)=0,         w_star(nx+2,ny+2,nz+1)=0, dp(nxp,nyp,nzp)=0
    real(8) :: rhs_x(nx+1,ny+2,nz+2)=0,          rhs_y(nx+2,ny+1,nz+2)=0,          rhs_z(nx+2,ny+2,nz+1)=0,   RHS_poisson(nxp,nyp,nzp)=0, RHS_poisson0(nxp*nyp*nzp)=0
    real(8) :: rhs_x_previous(nx+1,ny+2,nz+2)=0, rhs_y_previous(nx+2,ny+1,nz+2)=0, rhs_z_previous(nx+2,ny+2,nz+1)=0, dp_vec(nxp*nyp*nzp)=0, dp_lu(nxp,nyp,nzp)=0
    !real(8) :: rhs_x_previous0(nx+1,ny+2,nz+2)=0,rhs_y_previous0(nx+2,ny+1,nz+2)=0,rhs_z_previous0(nx+2,ny+2,nz+1)=0
    real(8) :: dpdx(nx+1,ny+2,nz+2)=0,           dpdy(nx+2,ny+1,nz+2)=0,           dpdz(nx+2,ny+2,nz+1)=0
    !real(8) :: f_term_x(nx+1,ny+2,nz+2)=0,       f_term_y(nx+2,ny+1,nz+2)=0,       f_term_z(nx+2,ny+2,nz+1)=0
    real(8) :: f_term_x=0,       f_term_y=0,       f_term_z=0

    real(8) :: bx_u_1(ny+2,nz+2)=0, bx_u_nx(ny+2,nz+2)=0, by_u_1(nx+1,nz+2)=0, by_u_ny(nx+1,nz+2)=0, bz_u_1(nx+1,ny+2)=0, bz_u_nz(nx+1,ny+2)=0
    real(8) :: bx_v_1(ny+1,nz+2)=0, bx_v_nx(ny+1,nz+2)=0, by_v_1(nx+2,nz+2)=0, by_v_ny(nx+2,nz+2)=0, bz_v_1(nx+2,ny+1)=0, bz_v_nz(nx+2,ny+1)=0
    real(8) :: bx_w_1(ny+2,nz+1)=0, bx_w_nx(ny+2,nz+1)=0, by_w_1(nx+2,nz+1)=0, by_w_ny(nx+2,nz+1)=0, bz_w_1(nx+2,ny+2)=0, bz_w_nz(nx+2,ny+2)=0
    real(8) :: bx_p_1(nyp,nzp)=0,   bx_p_nx(nyp,nzp)=0,   by_p_1(nxp,nzp)=0,   by_p_ny(nxp,nzp)=0,   bz_p_1(nxp,nyp)=0,   bz_p_nz(nxp,nyp)=0

    integer :: t_step, plot_step=20, slice=nz/2
    real(8) :: tGet
    integer, dimension (:), allocatable :: x_range0, y_range0, z_range0
    integer :: x_range1(nx)=[(i, i=2, nx+1)], y_range1(ny)=[(i, i=2, ny+1)], z_range1(nz)=[(i, i=2, nz+1)]
    type(csr), allocatable :: LHS_poisson
    type(coo), allocatable :: LHS_poisson_coo
    real(8), dimension (:,:), allocatable :: LHS_poisson_den

    real(8) :: temp01, temp02, temp03, temp04, temp05, temp06
    real(8), dimension (:), allocatable :: temp11, temp12, temp13, temp14, temp15, temp16
    real(8), dimension (:,:), allocatable :: temp21, temp22, temp23, temp24, temp25, temp26
    real(8), dimension (:,:,:), allocatable :: temp31, temp32, temp33, temp34, temp35, temp36
    real(4) :: ini_vel(3,256,256,256), ini_pr(256,256,256)
    real(8) :: div(nx,ny,nz)
    !real(8), dimension (:,:,:), allocatable :: vx, vy, vz, wx, wy, wz
    !real(8), dimension (:,:,:), allocatable :: conv_x, conv_y, conv_z, diff_x, diff_y, diff_z
    real(8) :: ux(nx+1,ny+2,nz+2), uy(nx+1,ny+1,nz+2), uz(nx+1,ny+2,nz+1)
    real(8) :: vx(nx+1,ny+1,nz+2), vy(nx+2,ny+1,nz+2), vz(nx+2,ny+1,nz+1)
    real(8) :: wx(nx+1,ny+2,nz+1), wy(nx+2,ny+1,nz+1), wz(nx+2,ny+2,nz+1)
    real(8) :: conv_x(nx,ny,nz), conv_y(nx,ny,nz), conv_z(nx,ny,nz)
    real(8) :: diff_x(nx,ny,nz), diff_y(nx,ny,nz), diff_z(nx,ny,nz)
    real(8) :: RHS_poisson_internal(nx,ny,nz)

    !INTEL mkl_pardiso
    type(MKL_PARDISO_HANDLE) pt(64)
    integer :: maxfct=1, mnum=1, mtype=11, phase=13, n=nxp*nyp*nzp, idum(nxp*nyp*nzp), nrhs=1, iparm(64)=0, msglvl=0, error=0

    !INTEL mkl_dft
    integer :: cstrides(4)=0, rstrides(4)=0
    type(DFTI_DESCRIPTOR), POINTER :: hand_f, hand_b
    real(8) :: dft_out_r(nx,ny,nz), poisson_eigv(nx,ny,nz)=0!, dft_in_r(nx,ny,nz)=0
    complex*16 :: dft_in_c(nx,ny,nz), dft_out_c(nx,ny,nz)
    integer :: status

    !post-processing
    type(gpf):: gp
    character(len=1000) :: string_var
    INTEGER :: funit, io_stat

    !system_clock
    REAL(8) :: system_clock_rate
    INTEGER :: c01,c02,c1,c2,cr,cm

    call OMP_set_dynamic(.true.)
    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)

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

    if (LU_poisson) then
        LHS_poisson=Poisson_LHS_staggered(nxp, nyp, nzp, dx2, dy2, dz2, pbc_x, pbc_y, pbc_z, dx, dy, dz)
        DO i = 1, 64
            pt(i)%DUMMY = 0
        END DO
        phase=12
        iparm(1) = 1 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(5) = 2 ! no user fill-in reducing permutation
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
        iparm(19) = -1 ! Output: Mflops for LU factorization

        print *, "**************************************"
        print *, "LU of LHS_poisson start..."
        CALL SYSTEM_CLOCK(c1)
        call pardiso (pt, maxfct, mnum, mtype, phase, n, LHS_poisson%value, LHS_poisson%ia, LHS_poisson%ja, &
            idum, nrhs, iparm, msglvl, RHS_poisson0, dp_vec, error)
        CALL SYSTEM_CLOCK(c2)
        print '(" LU of LHS_poisson completed: " (f6.4) " second")', (c2-c1)/system_clock_rate
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
    !status = DftiCreateDescriptor(hand_f, DFTI_DOUBLE, DFTI_REAL, 3, [nx,ny,nz])
    !status = DftiSetValue(hand_f, DFTI_BACKWARD_SCALE, 1.0d0/(nx*ny*nz))
    !status = DftiSetValue(hand_f, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    !status = DftiSetValue(hand_f, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
    !status = DftiSetValue(hand_f, DFTI_INPUT_STRIDES, rstrides)
    !status = DftiSetValue(hand_f, DFTI_OUTPUT_STRIDES, cstrides)
    status = DftiCreateDescriptor(hand_f, DFTI_DOUBLE, DFTI_COMPLEX, 3, [nx,ny,nz])
    status = DftiCommitDescriptor(hand_f)

    print *,"Configure DFTI descriptor for backward transform"
    !status = DftiCreateDescriptor(hand_b, DFTI_DOUBLE, DFTI_REAL, 3, [nx,ny,nz])
    !status = DftiSetValue(hand_b, DFTI_BACKWARD_SCALE, 1.0d0/(nx*ny*nz))
    !status = DftiSetValue(hand_b, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    !status = DftiSetValue(hand_b, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
    !status = DftiSetValue(hand_b, DFTI_INPUT_STRIDES, cstrides)
    !status = DftiSetValue(hand_b, DFTI_OUTPUT_STRIDES, rstrides)
    status = DftiCreateDescriptor(hand_b, DFTI_DOUBLE, DFTI_COMPLEX, 3, [nx,ny,nz])
    status = DftiSetValue(hand_b, DFTI_BACKWARD_SCALE, 1.0d0/(nx*ny*nz))
    status = DftiCommitDescriptor(hand_b)

    do i=1,nx
        poisson_eigv(i,:,:)=(sin(pi*(i-1)/nx)/dx)**2
    end do

    do i=1,ny
        poisson_eigv(:,i,:)=poisson_eigv(:,i,:)+(sin(pi*(i-1)/ny)/dy)**2
    end do

    do i=1,nz
        poisson_eigv(:,:,i)=poisson_eigv(:,:,i)+(sin(pi*(i-1)/nz)/dz)**2
    end do
    poisson_eigv(1,1,1)=1
    poisson_eigv=-poisson_eigv*4.0d0

    do t_step=0,time_length
        tGet=time_array(t_step)
        print *,''
        print '(" t_step " (i6))', t_step
        CALL SYSTEM_CLOCK(c01)

        if (t_step==0) then
            !load(['D:\Documents\HIT_32^3_data_staggered\HIT_32^3_init.mat'],'u','v','w','p');
            open(20, file="D:\Downloads\iso1024_256_cOrder.bin", form='unformatted', status='old', &
                access='direct', recl=3*256*256*256)
            READ(20, rec=1) ini_vel
            close(20)

            open(20, file="D:\Downloads\iso1024_256_cOrder_Pressure.bin", form='unformatted', status='old', &
                access='direct', recl=256*256*256)
            READ(20, rec=1) ini_pr
            close(20)

            temp01=256/nx; temp02=256/ny; temp03=256/nz;
            u(1:nx,2:ny+1,2:nz+1)=ini_vel(1,1:256:temp01,2:256:temp02,2:256:temp03)
            u(nx+1,:,:)=u(1,:,:)
            u(:,1,:)=u(:,ny+1,:); u(:,ny+2,:)=u(:,2,:)
            u(:,:,1)=u(:,:,nz+1); u(:,:,nz+2)=u(:,:,2)

            v(2:nx+1,1:ny,2:nz+1)=ini_vel(2,2:256:temp01,1:256:temp02,2:256:temp03)
            v(1,:,:)=v(nx+1,:,:); v(ny+2,:,:)=v(2,:,:)
            v(:,nx+1,:)=v(:,1,:);
            v(:,:,1)=v(:,:,nz+1); v(:,:,nz+2)=v(:,:,2)

            w(2:nx+1,2:ny+1,1:nz)=ini_vel(3,2:256:temp01,2:256:temp02,1:256:temp03)
            w(1,:,:)=w(nx+1,:,:); w(ny+2,:,:)=w(2,:,:)
            w(:,1,:)=w(:,ny+1,:); w(:,ny+2,:)=w(:,2,:)
            w(:,:,nz+1)=w(:,:,1);

            p(2:nx+1,2:ny+1,2:nz+1)=ini_pr(2:256:temp01,2:256:temp02,2:256:temp03)
            p(1,:,:)=p(nx+1,:,:); p(ny+2,:,:)=p(2,:,:)
            p(:,1,:)=p(:,ny+1,:); p(:,ny+2,:)=p(:,2,:)
            p(:,:,1)=p(:,:,nz+1); p(:,:,nz+2)=p(:,:,2)

            !call TGV(xu, yu, zu, 0.0d0, nu, u)
            !call TGV(xv, yv, zv, 0.0d0, nu, v=v)
            !call TGV(xp, yp, zp, 0.0d0, nu, p=p)
        else
            !!! Convection
            ux=avg(u,1,bc_x); uy=avg(u,2,0);    uz=avg(u,3,0);
            vx=avg(v,1,0);    vy=avg(v,2,bc_y); vz=avg(v,3,0);
            wx=avg(w,1,0);    wy=avg(w,2,0);    wz=avg(w,3,bc_z);

            CALL SYSTEM_CLOCK(c1)
            if (bc_x==1) then
                conv_x=diff(ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1)*ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1),1,1)/dx
                conv_x=conv_x + diff(uy(1:ubound(uy,1)-1,:,2:ubound(uy,3)-1)*vx(1:ubound(vx,1)-1,:,2:ubound(vx,3)-1),1,2)/dy
                conv_x=conv_x + diff(uz(1:ubound(uz,1)-1,2:ubound(uz,2)-1,:)*wx(1:ubound(wx,1)-1,2:ubound(wx,2)-1,:),1,3)/dz
                !conv_x=diff(ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1)*ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1),1,1)/dx + &
                !    diff(uy(1:ubound(uy,1)-1,:,2:ubound(uy,3)-1)*vx(1:ubound(vx,1)-1,:,2:ubound(vx,3)-1),1,2)/dy + &
                !    diff(uz(1:ubound(uz,1)-1,2:ubound(uz,2)-1,:)*wx(1:ubound(wx,1)-1,2:ubound(wx,2)-1,:),1,3)/dz
            else
                conv_x=diff(ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1)*ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1),1,1)/dx
                conv_x=conv_x + diff(uy(2:ubound(uy,1)-1,:,2:ubound(uy,3)-1)*vx(2:ubound(vx,1)-1,:,2:ubound(vx,3)-1),1,2)/dy
                conv_x=conv_x + diff(uz(2:ubound(uz,1)-1,2:ubound(uz,2)-1,:)*wx(2:ubound(wx,1)-1,2:ubound(wx,2)-1,:),1,3)/dz
                !conv_x=diff(ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1)*ux(:,2:ubound(ux,2)-1,2:ubound(ux,3)-1),1,1)/dx + &
                !    diff(uy(2:ubound(uy,1)-1,:,2:ubound(uy,3)-1)*vx(2:ubound(vx,1)-1,:,2:ubound(vx,3)-1),1,2)/dy + &
                !    diff(uz(2:ubound(uz,1)-1,2:ubound(uz,2)-1,:)*wx(2:ubound(wx,1)-1,2:ubound(wx,2)-1,:),1,3)/dz
            end if

            if (bc_y==1) then
                conv_y=diff(vx(:,1:ubound(vx,2)-1,2:ubound(vx,3)-1)*uy(:,1:ubound(uy,2)-1,2:ubound(uy,3)-1),1,1)/dx
                conv_y=conv_y + diff(vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1)*vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1),1,2)/dy
                conv_y=conv_y + diff(vz(2:ubound(vz,1)-1,1:ubound(vz,2)-1,:)*wy(2:ubound(wy,1)-1,1:ubound(wy,2)-1,:),1,3)/dz
                !conv_y=diff(vx(:,1:ubound(vx,2)-1,2:ubound(vx,3)-1)*uy(:,1:ubound(uy,2)-1,2:ubound(uy,3)-1),1,1)/dx + &
                !    diff(vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1)*vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1),1,2)/dy + &
                !    diff(vz(2:ubound(vz,1)-1,1:ubound(vz,2)-1,:)*wy(2:ubound(wy,1)-1,1:ubound(wy,2)-1,:),1,3)/dz
            else
                conv_y=diff(vx(:,2:ubound(vx,2)-1,2:ubound(vx,3)-1)*uy(:,2:ubound(uy,2)-1,2:ubound(uy,3)-1),1,1)/dx
                conv_y=conv_y + diff(vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1)*vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1),1,2)/dy
                conv_y=conv_y + diff(vz(2:ubound(vz,1)-1,2:ubound(vz,2)-1,:)*wy(2:ubound(wy,1)-1,2:ubound(wy,2)-1,:),1,3)/dz
                !conv_y=diff(vx(:,2:ubound(vx,2)-1,2:ubound(vx,3)-1)*uy(:,2:ubound(uy,2)-1,2:ubound(uy,3)-1),1,1)/dx + &
                !    diff(vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1)*vy(2:ubound(vy,1)-1,:,2:ubound(vy,3)-1),1,2)/dy + &
                !    diff(vz(2:ubound(vz,1)-1,2:ubound(vz,2)-1,:)*wy(2:ubound(wy,1)-1,2:ubound(wy,2)-1,:),1,3)/dz
            end if

            if (bc_z==1) then
                conv_z=diff(wx(:,2:ubound(wx,2)-1,1:ubound(wx,3)-1)*uz(:,2:ubound(uz,2)-1,1:ubound(uz,3)-1),1,1)/dx
                conv_z=conv_z + diff(wy(2:ubound(wy,1)-1,:,1:ubound(wy,3)-1)*vz(2:ubound(vz,1)-1,:,1:ubound(vz,3)-1),1,2)/dy
                conv_z=conv_z + diff(wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:)*wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:),1,3)/dz
                !conv_z=diff(wx(:,2:ubound(wx,2)-1,1:ubound(wx,3)-1)*uz(:,2:ubound(uz,2)-1,1:ubound(uz,3)-1),1,1)/dx + &
                !    diff(wy(2:ubound(wy,1)-1,:,1:ubound(wy,3)-1)*vz(2:ubound(vz,1)-1,:,1:ubound(vz,3)-1),1,2)/dy + &
                !    diff(wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:)*wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:),1,3)/dz
            else
                conv_z=diff(wx(:,2:ubound(wx,2)-1,2:ubound(wx,3)-1)*uz(:,2:ubound(uz,2)-1,2:ubound(uz,3)-1),1,1)/dx
                conv_z=conv_z + diff(wy(2:ubound(wy,1)-1,:,2:ubound(wy,3)-1)*vz(2:ubound(vz,1)-1,:,2:ubound(vz,3)-1),1,2)/dy
                conv_z=conv_z + diff(wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:)*wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:),1,3)/dz
                !conv_z=diff(wx(:,2:ubound(wx,2)-1,2:ubound(wx,3)-1)*uz(:,2:ubound(uz,2)-1,2:ubound(uz,3)-1),1,1)/dx + &
                !    diff(wy(2:ubound(wy,1)-1,:,2:ubound(wy,3)-1)*vz(2:ubound(vz,1)-1,:,2:ubound(vz,3)-1),1,2)/dy + &
                !    diff(wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:)*wz(2:ubound(wz,1)-1,2:ubound(wz,2)-1,:),1,3)/dz
            end if

            !!! Diffusion
            if (bc_x==1) then
                diff_x=diff2(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,bc_x)/dx2
                diff_x=diff_x + diff2(u(1:ubound(u,1)-1,:,2:ubound(u,3)-1),2,0)/dy2
                diff_x=diff_x + diff2(u(1:ubound(u,1)-1,2:ubound(u,2)-1,:),3,0)/dz2
                !diff_x=diff2(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,bc_x)/dx2+diff2(u(1:ubound(u,1)-1,:,2:ubound(u,3)-1),2,0)/dy2+diff2(u(1:ubound(u,1)-1,2:ubound(u,2)-1,:),3,0)/dz2
            else
                diff_x=diff2(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,bc_x)/dx2
                diff_x=diff_x + diff2(u(2:ubound(u,1)-1,:,2:ubound(u,3)-1),2,0)/dy2
                diff_x=diff_x + diff2(u(2:ubound(u,1)-1,2:ubound(u,2)-1,:),3,0)/dz2
                !diff_x=diff2(u(:,2:ubound(u,2)-1,2:ubound(u,3)-1),1,bc_x)/dx2+diff2(u(2:ubound(u,1)-1,:,2:ubound(u,3)-1),2,0)/dy2+diff2(u(2:ubound(u,1)-1,2:ubound(u,2)-1,:),3,0)/dz2
            end if

            if (bc_y==1) then
                diff_y=diff2(v(:,1:ubound(v,2)-1,2:ubound(v,3)-1),1,0)/dx2
                diff_y=diff_y + diff2(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),2,bc_y)/dy2
                diff_y=diff_y + diff2(v(2:ubound(v,1)-1,1:ubound(v,2)-1,:),3,0)/dz2
                !diff_y=diff2(v(:,1:ubound(v,2)-1,2:ubound(v,3)-1),1,0)/dx2+diff2(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),2,bc_y)/dy2+diff2(v(2:ubound(v,1)-1,1:ubound(v,2)-1,:),3,0)/dz2
            else
                diff_y=diff2(v(:,2:ubound(v,2)-1,2:ubound(v,3)-1),1,0)/dx2
                diff_y=diff_y + diff2(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),2,bc_y)/dy2
                diff_y=diff_y + diff2(v(2:ubound(v,1)-1,2:ubound(v,2)-1,:),3,0)/dz2
                !diff_y=diff2(v(:,2:ubound(v,2)-1,2:ubound(v,3)-1),1,0)/dx2+diff2(v(2:ubound(v,1)-1,:,2:ubound(v,3)-1),2,bc_y)/dy2+diff2(v(2:ubound(v,1)-1,2:ubound(v,2)-1,:),3,0)/dz2
            end if

            if (bc_z==1) then
                diff_z=diff2(w(:,2:ubound(w,2)-1,1:ubound(w,3)-1),1,0)/dx2
                diff_z=diff_z +diff2(w(2:ubound(w,1)-1,:,1:ubound(w,3)-1),2,0)/dy2
                diff_z=diff_z +diff2(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),3,bc_z)/dz2
                !diff_z=diff2(w(:,2:ubound(w,2)-1,1:ubound(w,3)-1),1,0)/dx2+diff2(w(2:ubound(w,1)-1,:,1:ubound(w,3)-1),2,0)/dy2+diff2(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),3,bc_z)/dz2
            else
                diff_z=diff2(w(:,2:ubound(w,2)-1,2:ubound(w,3)-1),1,0)/dx2
                diff_z=diff_z + diff2(w(2:ubound(w,1)-1,:,2:ubound(w,3)-1),2,0)/dy2
                diff_z=diff_z + diff2(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),3,bc_z)/dz2
                !diff_z=diff2(w(:,2:ubound(w,2)-1,2:ubound(w,3)-1),1,0)/dx2+diff2(w(2:ubound(w,1)-1,:,2:ubound(w,3)-1),2,0)/dy2+diff2(w(2:ubound(w,1)-1,2:ubound(w,2)-1,:),3,bc_z)/dz2
            end if
            CALL SYSTEM_CLOCK(c2)
            print '("    time1: ", (f6.4), " second")', (c2-c1)/system_clock_rate

            !!! Time-advancement
            f_term_x=0; f_term_y=0; f_term_z=0;

            CALL SYSTEM_CLOCK(c1)
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
            else if (timescheme=="AB2") then
                ! diff terms + conv terms
                rhs_x(x_range0,y_range1,z_range1)=nu*diff_x-conv_x;
                rhs_y(x_range1,y_range0,z_range1)=nu*diff_y-conv_y;
                rhs_z(x_range1,y_range1,z_range0)=nu*diff_z-conv_z;

                ! prediction
                if (t_step==1 .and. init) then
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
            end if
            CALL SYSTEM_CLOCK(c2)
            print '("    time2: ", (f6.4), " second")', (c2-c1)/system_clock_rate


            !print *, rhs_x_previous(4,4,2)
            !print *, rhs_y_previous(7,10,15)
            !print *, rhs_z_previous(3,4,2)

            call vel_bc_staggered(u_star,v_star,w_star,&
                bx_u_1,bx_u_nx,by_u_1,by_u_ny,bz_u_1,bz_u_nz,&
                bx_v_1,bx_v_nx,by_v_1,by_v_ny,bz_v_1,bz_v_nz,&
                bx_w_1,bx_w_nx,by_w_1,by_w_ny,bz_w_1,bz_w_nz,&
                bc_x,bc_y,bc_z);

            !print *, u_star(4,4,2)
            !print *, v_star(7,10,15)
            !print *, w_star(3,4,2)

            !RHS_poisson(2:ubound(RHS_poisson,1)-1,2:ubound(RHS_poisson,2)-1,2:ubound(RHS_poisson,3)-1) = &
            !    (diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx + &
            !    diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy + &
            !    diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz) / (dt0);
            RHS_poisson_internal = diff(u_star(:,2:ubound(u_star,2)-1,2:ubound(u_star,3)-1),1,1)/dx
            RHS_poisson_internal = RHS_poisson_internal + diff(v_star(2:ubound(v_star,1)-1,:,2:ubound(v_star,3)-1),1,2)/dy
            RHS_poisson_internal = RHS_poisson_internal + diff(w_star(2:ubound(w_star,1)-1,2:ubound(w_star,2)-1,:),1,3)/dz
            RHS_poisson_internal = RHS_poisson_internal/dt0

            if (LU_poisson) then
                RHS_poisson(2:ubound(RHS_poisson,1)-1,2:ubound(RHS_poisson,2)-1,2:ubound(RHS_poisson,3)-1)=RHS_poisson_internal
                if (pbc_x==1) then
                    RHS_poisson(1,:,:)=0;      RHS_poisson(ubound(RHS_poisson,1),:,:)=0;
                else if (pbc_x==2 .or. pbc_x==4) then
                    RHS_poisson(1,:,:)=bx_p_1; RHS_poisson(ubound(RHS_poisson,1),:,:)=bx_p_nx;
                else if (pbc_x==3) then
                    RHS_poisson(1,:,:)=bx_p_1; RHS_poisson(ubound(RHS_poisson,1),:,:)=bx_p_nx;
                end if
                if (pbc_y==1) then
                    RHS_poisson(:,1,:)=0;      RHS_poisson(:,ubound(RHS_poisson,2),:)=0;
                else if (pbc_y==2 .or. pbc_y==4) then
                    RHS_poisson(:,1,:)=by_p_1; RHS_poisson(:,ubound(RHS_poisson,2),:)=by_p_ny;
                else if (pbc_y==3) then
                    RHS_poisson(:,1,:)=by_p_1; RHS_poisson(:,ubound(RHS_poisson,2),:)=by_p_ny;
                end if
                if (pbc_z==1) then
                    RHS_poisson(:,:,1)=0;      RHS_poisson(:,:,ubound(RHS_poisson,3))=0;
                else if (pbc_y==2 .or. pbc_y==4) then
                    RHS_poisson(:,:,1)=bz_p_1; RHS_poisson(:,:,ubound(RHS_poisson,3))=bz_p_nz;
                else if (pbc_z==3) then
                    RHS_poisson(:,:,1)=bz_p_1; RHS_poisson(:,:,ubound(RHS_poisson,3))=bz_p_nz;
                end if

                !LHS_poisson_coo=from_csr(LHS_poisson)
                !LHS_poisson_den=LHS_poisson_coo%to_den()
                RHS_poisson0=[RHS_poisson]

                !DO i = 1, 64
                !    pt(i)%DUMMY = 0
                !END DO
                phase=33
                !iparm(1) = 1 ! no solver default
                !iparm(2) = 2 ! fill-in reordering from METIS
                !iparm(5) = 2 ! no user fill-in reducing permutation
                !iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
                !iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
                !iparm(19) = -1 ! Output: Mflops for LU factorization


                print *, "**************************************"
                print *, "   Solve Poisson (LU decomp) start..."
                CALL SYSTEM_CLOCK(c1)
                call pardiso (pt, maxfct, mnum, mtype, phase, n, LHS_poisson%value, LHS_poisson%ia, LHS_poisson%ja, &
                    idum, nrhs, iparm, msglvl, RHS_poisson0, dp_vec, error)
                CALL SYSTEM_CLOCK(c2)
                print '("    Solve Poisson (LU decomp) completed: " (f6.4) " second")', (c2-c1)/system_clock_rate
                print *, "**************************************"
                dp_lu=reshape(dp_vec,([nxp,nyp,nzp]))
                dp_lu=dp_lu-dp_lu(2,2,2)
            end if

            CALL SYSTEM_CLOCK(c1)
            !print *, "**************************************"
            !print *, "   Solve Poisson (FFT-based FD) start..."
            !dft_in_r=RHS_poisson(2:nxp-1,2:nyp-1,2:nzp-1)
            dft_out_c=cmplx(RHS_poisson_internal)

            !print *, "   Forward FFT..."
            !status = DftiComputeForward(hand_f, dft_in_r(:,1,1), dft_out_c(:,1,1))
            status = DftiComputeForward(hand_f, dft_out_c(:,1,1))

            dft_out_c=dft_out_c/poisson_eigv

            !print *, "   Backward FFT..."
            !status = DftiComputeBackward(hand_b, dft_out_c(:,1,1), dft_out_r(:,1,1))
            status = DftiComputeBackward(hand_b, dft_out_c(:,1,1))
            dft_out_r=real(dft_out_c)
            dft_out_r=dft_out_r-dft_out_r(1,1,1)
            dp(2:nxp-1,2:nyp-1,2:nzp-1)=dft_out_r
            dp(1,:,:)=dp(nxp-1,:,:); dp(nxp,:,:)=dp(2,:,:)
            dp(:,1,:)=dp(:,nxp-1,:); dp(:,nxp,:)=dp(:,2,:)
            dp(:,:,1)=dp(:,:,nxp-1); dp(:,:,nxp)=dp(:,:,2)

            CALL SYSTEM_CLOCK(c2)
            print '("    Solve Poisson (FFT-based FD) completed: ", (f6.4), " second")', (c2-c1)/system_clock_rate
            !print *, "**************************************"

            !print *, maxval(abs(dp_dft-dp))

            !OPEN(10, file="poisson_eq.dat", form="unformatted")
            !WRITE(10) [LHS_poisson_den]
            !WRITE(10) [RHS_poisson]
            !WRITE(10) [dp]
            !WRITE(10) [dp_dft]
            !CLOSE(10)
            !dp=LHS_poisson\RHS_poisson(:);
            !dp=reshape(dp,nxp,nyp,nzp);
            !dp=dp_lu

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
        print '("Complete: ", f6.4, " second. MAX Div: ", e13.6)', (c02-c01)/system_clock_rate, maxval(abs(div))

        if (mod(t_step,plot_step)==0 .and. .true.) then
            plot_block: block

                !real(8), dimension (:,:), allocatable :: xx, yy

                !call meshgrid( xx, yy, yu(2:ubound(yu,1)-1), xu )
                write (string_var,'("iso256_t_step_",i0.6)') t_step
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
                call gp%title("DIV")
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

    !open(20, file="D:\Downloads\iso1024_256_cOrder.bin", form='unformatted', status='old', &
    !    access='direct', recl=1)
    !mm=0
    !do i=1,11
    !    do j=1,13
    !        do k=1,12
    !            do ll=1,3
    !                mm=mm+1
    !                !READ(20, rec=mm) temp41(i,j,k,ll)
    !                READ(20, rec=mm) temp42(ll,k,j,i)
    !            end do
    !        end do
    !    end do
    !end do


    !CALL SYSTEM_CLOCK(c01)
    !do i=1,40
    !    dpdx=diff(dp,1,1)/dx
    !    dpdy=diff(dp,1,2)/dy
    !    dpdz=diff(dp,1,3)/dz
    !    u=-dt0*dpdx+u_star
    !    v=-dt0*dpdy+v_star
    !    w=-dt0*dpdz+w_star
    !end do
    !CALL SYSTEM_CLOCK(c02)
    !print '(" time 1: ", (f6.4), " second")', (c02-c01)/system_clock_rate

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !if (.false.) then
    !    test_block: block
    !        !real(8), dimension (:,:,:), allocatable :: temp1, temp2, temp3, temp4, temp5, temp6
    !
    !        real(8), dimension(:), allocatable :: a_coo, b_coo
    !        integer, dimension(:), allocatable :: a_rowind, a_colind, b_rowind, b_colind
    !        integer, dimension(:), allocatable :: a_shape, b_shape
    !        real(8), dimension(:,:), allocatable :: c_mat
    !        type(coo) :: coo1, coo2, coo3
    !        type(csr) :: csr1, csr2, csr3
    !        integer :: info=-111, job(8)=0
    !
    !        !INTEL mkl_dft
    !        type(DFTI_DESCRIPTOR), POINTER :: hand
    !        real(8) :: dft_in_r_test(10)=0, dft_out_r_test(10)
    !        complex*16 :: dft_in_c_test(10), dft_out_c_test(10)
    !        integer :: stat
    !
    !        !fftw3
    !        integer :: fwd, bwd
    !
    !        !real(8), pointer :: temp11(), temp12
    !
    !        CALL SYSTEM_CLOCK(c01)
    !        do i=1,10
    !            temp31=diff(p, 1, 1)
    !        end do
    !        CALL SYSTEM_CLOCK(c02)
    !        print '(" timing 1: ", (f6.2), " second")', (c02-c01)/system_clock_rate
    !
    !        CALL SYSTEM_CLOCK(c02)
    !        do i=1,10
    !            temp32=diff_old(p, 1, 1)
    !        end do
    !        CALL SYSTEM_CLOCK(c02)
    !        print '(" timing 2: ", (f6.2), " second")', (c02-c01)/system_clock_rate
    !
    !        temp33=diff(p, 1, 2)
    !        temp34=diff_old(p, 1, 2)
    !        temp35=diff(p, 1, 3)
    !        temp36=diff_old(p, 1, 3)
    !
    !        print *, 'diff'
    !        print *, all(temp31==temp32)
    !        print *, all(temp33==temp34)
    !        print *, all(temp35==temp36)
    !
    !        print *, 'diff2'
    !        temp31=diff2(u,1,1);
    !        temp32=diff2_old(u,1,1);
    !        print *, all(temp31==temp32)
    !        temp31=diff2(u,1,0);
    !        temp32=diff2_old(u,1,0);
    !        print *, all(temp31==temp32)
    !        temp31=diff2(u,2,1);
    !        temp32=diff2_old(u,2,1);
    !        print *, all(temp31==temp32)
    !        temp31=diff2(u,2,0);
    !        temp32=diff2_old(u,2,0);
    !        print *, all(temp31==temp32)
    !        temp31=diff2(u,3,1);
    !        temp32=diff2_old(u,3,1);
    !        print *, all(temp31==temp32)
    !        temp31=diff2(u,3,0);
    !        temp32=diff2_old(u,3,0);
    !        print *, all(temp31==temp32)
    !
    !        print *, 'avg'
    !        temp31=avg(u,1,1);
    !        temp32=avg_old(u,1,1);
    !        print *, all(temp31==temp32)
    !        temp31=avg(u,1,0);
    !        temp32=avg_old(u,1,0);
    !        print *, all(temp31==temp32)
    !        temp31=avg(u,2,1);
    !        temp32=avg_old(u,2,1);
    !        print *, all(temp31==temp32)
    !        temp31=avg(u,2,0);
    !        temp32=avg_old(u,2,0);
    !        print *, all(temp31==temp32)
    !        temp31=avg(u,3,1);
    !        temp32=avg_old(u,3,1);
    !        print *, all(temp31==temp32)
    !        temp31=avg(u,3,0);
    !        temp32=avg_old(u,3,0);
    !        print *, all(temp31==temp32)
    !
    !        !print *, temp1(4,4,2)
    !        !print *, temp2(1,2,3)
    !        !print *, temp3(1,2,3)
    !        !print *, temp4(1,2,3)
    !        !print *, temp5(2,4,1)
    !        !print *, temp6(1,1,5)
    !
    !        !call coo_reshape(a_coo, a_rowind, a_colind, a_shape, (/1,6/))
    !        a_coo=(/1,-1,-3,-2,0,4,6,4,-4,2,7,8,-5/); a_rowind=(/1,1,1,2,2,3,3,3,4,4,4,5,5/); a_colind=(/1,2,3,1,2,3,4,5,1,3,4,2,5/); a_shape=(/5,5/)
    !        b_coo=(/1,2,3,4/); b_rowind=(/1,1,2,2/); b_colind=(/1,2,1,2/); b_shape=(/2,2/)
    !        coo1=coo_init_or_clean(a_coo,a_rowind,a_colind,a_shape)
    !        coo2=coo_init_or_clean(b_coo,b_rowind,b_colind,b_shape)
    !        coo3=kron(coo1,coo2)
    !        c_mat=coo1%to_den()
    !        call print_mat(c_mat)
    !        c_mat=coo2%to_den()
    !        call print_mat(c_mat)
    !        c_mat=coo3%to_den()
    !        call print_mat(c_mat)
    !
    !        csr1=coo1%to_csr()
    !        coo2=from_csr(csr1)
    !
    !        !LHS_poisson=Poisson_LHS_staggered(nxp, nyp, nzp, 0.01, 0.01, 0.01, pbc_x, pbc_y, pbc_z, 0.1, 0.1, 0.1)
    !        !coo2=from_csr(LHS_poisson)
    !        !c_mat=coo2%to_den()
    !        !call print_spmat(c_mat)
    !
    !        a_rowind=[1,2,3,4,5,6,7,8]
    !        a_colind=[2,4,6]
    !        print *, (a_rowind==a_colind)
    !
    !        dft_in_r_test=RHS_poisson(1:10,4,4)
    !        dft_in_c_test=cmplx(dft_in_r_test)
    !        print *, "dft_in"
    !        print *, dft_in_c_test
    !        print *, "**************************************"
    !
    !        fwd=0
    !        call dfftw_plan_dft_1d(fwd, 10, dft_in_c_test, dft_out_c_test, FFTW_FORWARD, FFTW_ESTIMATE)
    !        CALL SYSTEM_CLOCK(c1)
    !        call dfftw_execute(fwd)
    !        CALL SYSTEM_CLOCK(c2)
    !        print *, "**************************************"
    !        print '("    dft_forward FFTW3: ", (f6.4), " second")', (c2-c1)/system_clock_rate
    !        print *, dft_out_c_test
    !        print *, "**************************************"
    !
    !        dft_in_r_test=RHS_poisson(1:10,4,4)
    !        dft_in_c_test=cmplx(dft_in_r_test)
    !        ! Configure a Descriptor
    !        hand => null()
    !        !3D MKL dft = fft_z( fft_y (fft_x) ), same in MATLAB
    !        status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 1, 10)
    !        status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, 1.0d0/(10.0d0))
    !        status = DftiCommitDescriptor(hand)
    !        CALL SYSTEM_CLOCK(c1)
    !        status = DftiComputeForward(hand, dft_in_c_test)
    !        CALL SYSTEM_CLOCK(c2)
    !        print *, "**************************************"
    !        print '("    dft_forward MKL: ", (f6.4), " second")', (c2-c1)/system_clock_rate
    !        print *, dft_in_c_test
    !        print *, "**************************************"
    !
    !        CALL SYSTEM_CLOCK(c1)
    !        status = DftiComputeBackward(hand, dft_in_c_test)
    !        CALL SYSTEM_CLOCK(c2)
    !        print *, "**************************************"
    !        print '("    dft_backward MKL: ", (f6.4), " second")', (c2-c1)/system_clock_rate
    !        print *, dft_in_c_test
    !        print *, "**************************************"
    !        print *, maxval(abs(dft_in_c_test-dft_in_r_test))
    !
    !    end block test_block
    !end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) "Press anykey to continue..."
    !read (*,*)

    end program NS_3D_Staggered_dp

