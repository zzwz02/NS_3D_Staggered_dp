    subroutine test_subroutine

    implicit none

    include 'fftw/fftw3.f'

    ! Variables
    real(8), parameter :: pi = 3.1415926535897932_8
    integer :: i=0,j=0,k=0,ll=0,mm=0

    logical, parameter :: init=.true.
    character *4 :: timescheme="AB2"
    !integer, parameter :: timescheme=2 ! 1 Euler; 2 AB2
    real(8), parameter :: Re=1.0d0, nu=0.0008d0, t_start=-1.0d0, t_end=20.0d0, dt0=0.005d0
    integer, parameter :: nx_file=256
    integer, parameter :: nx=16, ny=nx, nz=nx, nxp=nx+2, nyp=ny+2, nzp=nz+2, sub_tstep=1
    logical, parameter :: LU_poisson=(nxp*nyp*nzp<=34**3)
    real(8), parameter :: lx=2.0d0*pi, ly=2.0d0*pi, lz=2.0d0*pi, dx=lx/nx, dy=ly/ny, dz=lz/nz, dx2=dx*dx, dy2=dy*dy, dz2=dz*dz
    real(8), parameter :: xu(nx+1)=[(i*dx, i=0, nx)],          yu(ny+2)=[((i+0.5)*dy, i=-1, ny)],   zu(nz+2)=[((i+0.5)*dz, i=-1, nz)]
    real(8), parameter :: xv(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yv(ny+1)=[(i*dy, i=0, ny)],          zv(nz+2)=[((i+0.5d0)*dz, i=-1, nz)]
    real(8), parameter :: xw(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yw(ny+2)=[((i+0.5d0)*dy, i=-1, ny)], zw(nz+1)=[(i*dz, i=0, nz)]
    real(8), parameter :: xp(nx+2)=[((i+0.5d0)*dx, i=-1, nx)], yp(ny+2)=[((i+0.5d0)*dy, i=-1, ny)], zp(nz+2)=[((i+0.5d0)*dz, i=-1, nz)]

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall)
    ! pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    integer, parameter :: bc_x=1, bc_y=1, bc_z=1, pbc_x=1, pbc_y=1, pbc_z=1

    integer, parameter :: time_length=nint((t_end-t_start)/(dt0))
    real(8), dimension (0:time_length), parameter :: time_array=[(t_start+dt0*i, i=0, time_length)]
    real(8) :: u(nx+1,ny+2,nz+2)=0,              v(nx+2,ny+1,nz+2)=0,              w(nx+2,ny+2,nz+1)=0,       p(nxp,nyp,nzp)=0
    real(8) :: dpdx(nx+1,ny+2,nz+2)=0,           dp_vec(nxp*nyp*nzp)=0,            dp_lu(nxp,nyp,nzp)=0
    real(8) ::  RHS_poisson(nxp,nyp,nzp)=0,      RHS_poisson0(nxp*nyp*nzp)=0

    type(csr), allocatable :: LHS_poisson
    type(coo), allocatable :: LHS_poisson_coo
    real(8), dimension (:,:), allocatable :: LHS_poisson_den

    real(8), dimension(:), allocatable :: a_coo, b_coo
    integer, dimension(:), allocatable :: a_rowind, a_colind, b_rowind, b_colind
    integer, dimension(:), allocatable :: a_shape, b_shape
    real(8), dimension(:,:), allocatable :: c_mat
    type(coo) :: coo1, coo2, coo3
    type(csr) :: csr1, csr2, csr3
    integer :: info=-111, job(8)=0

    integer(8) :: sizeof_record, sizeof_record_sub, tempi1, tempi2
    real(8) :: temp01, temp02, temp03, temp04, temp05, temp06
    real(8), dimension (:), allocatable :: temp11, temp12, temp13, temp14, temp15, temp16
    real(8), dimension (:,:), allocatable :: temp21, temp22, temp23, temp24, temp25, temp26
    real(8), dimension (:,:,:), allocatable :: temp31, temp32, temp33, temp34, temp35, temp36
    real(4) :: ini_vel(3,nx_file,nx_file,nx_file), ini_pr(nx_file,nx_file,nx_file)
    real(8) :: div(nx,ny,nz)

    !INTEL mkl_pardiso
    type(MKL_PARDISO_HANDLE) pt(64)
    integer :: maxfct=1, mnum=1, mtype=11, phase=13, n=nxp*nyp*nzp, idum(nxp*nyp*nzp), nrhs=1, iparm(64)=0, msglvl=0, error=0

    !INTEL mkl_dft
    type(DFTI_DESCRIPTOR), POINTER :: hand
    real(8) :: dft_in_r_test(10)=0, dft_out_r_test(10)
    complex*16 :: dft_in_c_test(10), dft_out_c_test(10)
    integer :: stat

    ! INTEL mkl_tt
    integer, parameter :: nx1=1024, ny1=1024
    real(8), parameter :: L1=1.0d0, dx1=L1/nx1, dy1=L1/ny1, dx21=dx1*dx1, dy21=dy1*dy1
    real(8), parameter :: xp1(nx1)=[((i-0.5d0)*dx1, i=1, nx1)], yp1(ny1)=[((i-0.5d0)*dy1, i=1, ny1)]
    integer :: tt_type, stat1, ipar_x(128), ipar_y(128)
    real(8) :: dpar_x(5*nx1/2+2)=0, dpar_y(5*ny1/2+2)=0, rhs_tt(nx1+1, ny1+1)=0, eig_tt(nx1, ny1), rhs_tt2(nx1+1, ny1+1)=0, rhs_tt1(nx1+1, ny1+1)=0
    type(dfti_descriptor), pointer :: handle_tx, handle_ty

    !fftw3
    integer :: fwd, bwd

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
    print '(" LU of LHS_poisson completed: ", F8.4, " second")', (c2-c1)/system_clock_rate
    print *, "**************************************"

    call TGV(xu, yu, zu, 0.0d0, nu, u)
    call TGV(xv, yv, zv, 0.0d0, nu, v=v)
    call TGV(xp, yp, zp, 0.0d0, nu, p=p)

    CALL SYSTEM_CLOCK(c01)
    do i=1,40
        dpdx=diff(p,1,1)/dx
        u=-dt0*dpdx+u
    end do
    CALL SYSTEM_CLOCK(c02)
    print '(" time 1: ", F8.4, " second")', (c02-c01)/system_clock_rate


    CALL SYSTEM_CLOCK(c01)
    do i=1,10
        temp31=diff(p, 1, 1)
    end do
    CALL SYSTEM_CLOCK(c02)
    print '(" timing 1: ", F8.4, " second")', (c02-c01)/system_clock_rate

    CALL SYSTEM_CLOCK(c02)
    do i=1,10
        temp32=diff_old(p, 1, 1)
    end do
    CALL SYSTEM_CLOCK(c02)
    print '(" timing 2: ", F8.4, " second")', (c02-c01)/system_clock_rate

    temp33=diff(p, 1, 2)
    temp34=diff_old(p, 1, 2)
    temp35=diff(p, 1, 3)
    temp36=diff_old(p, 1, 3)

    print *, 'diff'
    print *, all(temp31==temp32)
    print *, all(temp33==temp34)
    print *, all(temp35==temp36)

    print *, 'diff2'
    temp31=diff2(u,1,1);
    temp32=diff2_old(u,1,1);
    print *, all(temp31==temp32)
    temp31=diff2(u,1,0);
    temp32=diff2_old(u,1,0);
    print *, all(temp31==temp32)
    temp31=diff2(u,2,1);
    temp32=diff2_old(u,2,1);
    print *, all(temp31==temp32)
    temp31=diff2(u,2,0);
    temp32=diff2_old(u,2,0);
    print *, all(temp31==temp32)
    temp31=diff2(u,3,1);
    temp32=diff2_old(u,3,1);
    print *, all(temp31==temp32)
    temp31=diff2(u,3,0);
    temp32=diff2_old(u,3,0);
    print *, all(temp31==temp32)

    print *, 'avg'
    temp31=avg(u,1,1);
    temp32=avg_old(u,1,1);
    print *, all(temp31==temp32)
    temp31=avg(u,1,0);
    temp32=avg_old(u,1,0);
    print *, all(temp31==temp32)
    temp31=avg(u,2,1);
    temp32=avg_old(u,2,1);
    print *, all(temp31==temp32)
    temp31=avg(u,2,0);
    temp32=avg_old(u,2,0);
    print *, all(temp31==temp32)
    temp31=avg(u,3,1);
    temp32=avg_old(u,3,1);
    print *, all(temp31==temp32)
    temp31=avg(u,3,0);
    temp32=avg_old(u,3,0);
    print *, all(temp31==temp32)

    a_coo=(/1,-1,-3,-2,0,4,6,4,-4,2,7,8,-5/); a_rowind=(/1,1,1,2,2,3,3,3,4,4,4,5,5/); a_colind=(/1,2,3,1,2,3,4,5,1,3,4,2,5/); a_shape=(/5,5/)
    b_coo=(/1,2,3,4/); b_rowind=(/1,1,2,2/); b_colind=(/1,2,1,2/); b_shape=(/2,2/)
    coo1=coo_init_or_clean(a_coo,a_rowind,a_colind,a_shape)
    coo2=coo_init_or_clean(b_coo,b_rowind,b_colind,b_shape)
    coo3=kron(coo1,coo2)
    c_mat=coo1%to_den()
    call print_mat(c_mat)
    c_mat=coo2%to_den()
    call print_mat(c_mat)
    c_mat=coo3%to_den()
    call print_mat(c_mat)

    csr1=coo1%to_csr()
    coo2=from_csr(csr1)

    LHS_poisson=Poisson_LHS_staggered(nxp, nyp, nzp, 0.01d0, 0.01d0, 0.01d0, pbc_x, pbc_y, pbc_z, 0.1d0, 0.1d0, 0.1d0)
    coo2=from_csr(LHS_poisson)
    c_mat=coo2%to_den()
    !call print_spmat(c_mat)

    a_rowind=[1,2,3,4,5,6,7,8]
    a_colind=[2,4,6]
    print *, (a_rowind==a_colind)

    dft_in_r_test=RHS_poisson(1:10,4,4)
    dft_in_c_test=dcmplx(dft_in_r_test)
    print *, "dft_in"
    print *, dft_in_c_test
    print *, "**************************************"

    fwd=0
    call dfftw_plan_dft_1d(fwd, 10, dft_in_c_test, dft_out_c_test, FFTW_FORWARD, FFTW_ESTIMATE)
    CALL SYSTEM_CLOCK(c1)
    call dfftw_execute(fwd)
    CALL SYSTEM_CLOCK(c2)
    print *, "**************************************"
    print '("    dft_forward FFTW3: ", F8.4, " second")', (c2-c1)/system_clock_rate
    print *, dft_out_c_test
    print *, "**************************************"

    dft_in_r_test=RHS_poisson(1:10,4,4)
    dft_in_c_test=dcmplx(dft_in_r_test)
    ! Configure a Descriptor
    hand => null()
    !3D MKL dft = fft_z( fft_y (fft_x) ), same in MATLAB
    stat = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 1, 10)
    stat = DftiSetValue(hand, DFTI_BACKWARD_SCALE, 1.0d0/(10.0d0))
    stat = DftiCommitDescriptor(hand)
    CALL SYSTEM_CLOCK(c1)
    stat = DftiComputeForward(hand, dft_in_c_test)
    CALL SYSTEM_CLOCK(c2)
    print *, "**************************************"
    print '("    dft_forward MKL: ", F8.4, " second")', (c2-c1)/system_clock_rate
    print *, dft_in_c_test
    print *, "**************************************"

    CALL SYSTEM_CLOCK(c1)
    stat = DftiComputeBackward(hand, dft_in_c_test)
    CALL SYSTEM_CLOCK(c2)
    print *, "**************************************"
    print '("    dft_backward MKL: ", F8.4, " second")', (c2-c1)/system_clock_rate
    print *, dft_in_c_test
    print *, "**************************************"
    print *, maxval(abs(dft_in_c_test-dft_in_r_test))

    print *, "**************************************"
    print *,"D-D:   DST-1, Forward sine transform"
    print *,"D-N:   DST-3, Forward staggered sine transform"
    print *,"Ds-Ns: DST-4, Forward staggered2 sine transform"
    print *,"N-N:   DCT-1, Forward cosine transform"
    print *,"N-D:   DCT-3, Forward staggered cosine transform"
    print *,"Ns-Ds: DCT-4, Forward staggered2 cosine transform"
    print *, "**************************************"

    print *,"Example:"
    print *,"  x: Ds-Ds b.c. (DST-2)"
    print *,"    DST-2 == flipping sign of every other input, DCT-2 (or iDCT-3 with 2/N factor), reversing the order of output"
    print *,"    IDST-2 == reversing the order of input, iDCT-2 (or DCT-3 with 2/N factor), flipping sign of every other output"
    print *,"  y: Ns-Ns b.c. (DCT-2)"
    print *,"    DCT-2 == iDCT-3 with 2/N factor"
    print *,"    IDCT-2 == DCT-3"

    call d_init_trig_transform(nx1, MKL_STAGGERED_COSINE_TRANSFORM, ipar_x, dpar_x, stat)
    call d_init_trig_transform(ny1, MKL_STAGGERED_COSINE_TRANSFORM, ipar_y, dpar_y, stat)
    call d_commit_trig_transform(rhs_tt(:,1),handle_tx,ipar_x,dpar_x,stat)
    call d_commit_trig_transform(rhs_tt(1,:),handle_ty,ipar_y,dpar_y,stat)
    
    do j=1,ny1
        rhs_tt(1:nx1,j)=-2*(2*yp1(j)**3-3*yp1(j)**2+1)+6*(1-xp1**2)*(2*yp1(j)-1);
    end do
    
    do i=1,nx1
        eig_tt(i,:)=(sin(i*pi/2/nx1))**2/dx21
    end do

    do j=1,ny1
        eig_tt(:,j)=eig_tt(:,j)+(sin((j-1)*pi/2/ny1))**2/dy21
    end do
    eig_tt=-4.0d0*eig_tt

    !rhs_tt(1:nx)=[(i*sin(dble(i)), i=1,nx)]
    rhs_tt(1,1:ny1)  =rhs_tt(1,1:ny1)  -2*(2*yp1**3-3*yp1**2+1)/dx21 !b.c condition
    rhs_tt(nx1,1:ny1)=rhs_tt(nx1,1:ny1)-2*(0)/dx21 !b.c condition
    rhs_tt(1:nx1,1)  =rhs_tt(1:nx1,1)  +(1)/dy1 !b.c condition
    rhs_tt(1:nx1,ny1)=rhs_tt(1:nx1,ny1)-(0)/dy1 !b.c condition
    rhs_tt1=rhs_tt
    rhs_tt2=rhs_tt
    
    CALL SYSTEM_CLOCK(c1)
    do i=1,nx1
        call d_backward_trig_transform(rhs_tt(i,:),handle_ty,ipar_y,dpar_y,stat)
    end do
    !rhs_tt(:,1)=rhs_tt(:,1)/sqrt(dble(ny1))
    !rhs_tt(:,2:ny1)=rhs_tt(:,2:ny1)/sqrt(ny1/2.0d0)

    do j=1,ny1
        rhs_tt(1:nx1:2,j)=-rhs_tt(1:nx1:2,j)
        call d_backward_trig_transform(rhs_tt(:,j),handle_tx,ipar_x,dpar_x,stat)
        rhs_tt(1:nx1,j)=rhs_tt(nx1:1:-1,j)
    end do
    !rhs_tt(:,:)=rhs_tt(:,:)/sqrt(nx1/2.0d0)
    
    rhs_tt(1:nx1,1:ny1)=rhs_tt(1:nx1,1:ny1)/eig_tt!(nx1:1:-1,:)
    
    do j=1,ny1
        rhs_tt(1:nx1,j)=rhs_tt(nx1:1:-1,j)
        call d_forward_trig_transform(rhs_tt(:,j),handle_tx,ipar_x,dpar_x,stat)
        rhs_tt(1:nx1:2,j)=-rhs_tt(1:nx1:2,j)
    end do
    !rhs_tt(:,:)=rhs_tt(:,:)*sqrt(nx1/2.0d0)
    
    do i=1,nx1
        call d_forward_trig_transform(rhs_tt(i,:),handle_ty,ipar_y,dpar_y,stat)
    end do
    !rhs_tt(:,2:ny1)=rhs_tt(:,2:ny1)*sqrt(ny1/2.0d0)
    !rhs_tt(:,1)=rhs_tt(:,1)*sqrt(dble(ny1))
    CALL SYSTEM_CLOCK(c2)
    print *, "**************************************"
    print '("    dct_2D Poisson MKL1: ", F8.4, " second")', (c2-c1)/system_clock_rate
    
    CALL SYSTEM_CLOCK(c1)
    do i=1,nx1
        call d_backward_trig_transform(rhs_tt1(i,:),handle_ty,ipar_y,dpar_y,stat)
    end do
    !rhs_tt1(:,1)=rhs_tt1(:,1)/sqrt(dble(ny1))
    !rhs_tt1(:,2:ny1)=rhs_tt1(:,2:ny1)/sqrt(ny1/2.0d0)

    rhs_tt1(1:nx1:2,:)=-rhs_tt1(1:nx1:2,:)
    do j=1,ny1
        call d_backward_trig_transform(rhs_tt1(:,j),handle_tx,ipar_x,dpar_x,stat)
    end do
    rhs_tt1(1:nx1,:)=rhs_tt1(nx1:1:-1,:)
    !rhs_tt1(:,:)=rhs_tt1(:,:)/sqrt(nx1/2.0d0)
    
    rhs_tt1(1:nx1,1:ny1)=rhs_tt1(1:nx1,1:ny1)/eig_tt
    
    rhs_tt1(1:nx1,:)=rhs_tt1(nx1:1:-1,:)
    do j=1,ny1
        call d_forward_trig_transform(rhs_tt1(:,j),handle_tx,ipar_x,dpar_x,stat)
    end do
    rhs_tt1(1:nx1:2,:)=-rhs_tt1(1:nx1:2,:)
    !rhs_tt1(:,:)=rhs_tt1(:,:)*sqrt(nx1/2.0d0)
    
    do i=1,nx1
        call d_forward_trig_transform(rhs_tt1(i,:),handle_ty,ipar_y,dpar_y,stat)
    end do
    !rhs_tt1(:,2:ny1)=rhs_tt1(:,2:ny1)*sqrt(ny1/2.0d0)
    !rhs_tt1(:,1)=rhs_tt1(:,1)*sqrt(dble(ny1))
    CALL SYSTEM_CLOCK(c2)
    print *, "**************************************"
    print '("    dct_2D Poisson MKL2: ", F8.4, " second")', (c2-c1)/system_clock_rate
    
    CALL SYSTEM_CLOCK(c1)
    do i=1,nx1
        call DST_DCT_2(rhs_tt2(i,:), ny1, handle_ty, ipar_y, dpar_y, DCT2=.true., forward=.true.)
    end do

    do j=1,ny1
        call DST_DCT_2(rhs_tt2(:,j), nx1, handle_tx, ipar_x, dpar_x, DST2=.true., forward=.true.)
    end do
    
    rhs_tt2(1:nx1,1:ny1)=rhs_tt2(1:nx1,1:ny1)/eig_tt
    
    do j=1,ny1
        call DST_DCT_2(rhs_tt2(:,j), nx1, handle_tx, ipar_x, dpar_x, DST2=.true., backward=.true.)
    end do
    
    do i=1,nx1
        call DST_DCT_2(rhs_tt2(i,:), ny1, handle_ty, ipar_y, dpar_y, DCT2=.true., backward=.true.)
    end do
    CALL SYSTEM_CLOCK(c2)
    print *, "**************************************"
    print '("    dct_2D Poisson MKL3: ", F8.4, " second")', (c2-c1)/system_clock_rate
    
    !rhs_tt(2:nx+1:2)=-rhs_tt(2:nx+1:2)
    !eig=[(-4/dx2*(sin(i*pi/2/nx))**2, i=1,nx)]
    !call d_init_trig_transform(nx, MKL_STAGGERED_COSINE_TRANSFORM, ipar, dpar, stat)
    !call d_commit_trig_transform(rhs_tt,handle_tt,ipar,dpar,stat)
    !
    !CALL D_BACKWARD_TRIG_TRANSFORM(rhs_tt,handle_tt,ipar,dpar,stat)
    !rhs_tt(1:nx)=rhs_tt(1:nx)/eig(nx:1:-1)
    !CALL D_FORWARD_TRIG_TRANSFORM(rhs_tt,handle_tt,ipar,dpar,stat)
    !rhs_tt(2:nx+1:2)=-rhs_tt(2:nx+1:2)


    end subroutine test_subroutine