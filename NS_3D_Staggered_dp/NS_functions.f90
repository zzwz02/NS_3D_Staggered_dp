    module NS_functions

    use coo_mod
    use other_utility

    implicit none

    include "mkl_spblas.fi"

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine TGV(x, y, z, t, nu, u, v, w, p)
    implicit none

    real(8), value :: t, nu
    real(8), dimension (:), intent(in) :: x, y, z
    real(8), dimension (:,:,:), intent(inout), optional :: u, v, w, p

    integer :: nx, ny, nz, i, j, k
    
    !if (size(A,1)==1) A=transpose(A)

    nx=size(x)
    ny=size(y)
    nz=size(z)

    nu=exp(-2*nu*t)

    do j=1, ny
        do i=1, nx
            if (Present (u)) then
                u(i,j,:)=+cos(x(i))*sin(y(j))*nu
            end if
            if (Present (v)) then
                v(i,j,:)=-sin(x(i))*cos(y(j))*nu
            end if
            if (Present (p)) then
                p(i,j,:)=-0.25d0*(cos(2.0d0*x(i))+cos(2.0d0*y(j)))*nu*nu
            end if
        end do
    end do

    end subroutine TGV

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine vel_bc_staggered(u_star, v_star, w_star,&
        bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz,&
        bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz,&
        bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz,&
        bc_x, bc_y, bc_z)

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall)
    ! pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell

    implicit none
    integer, intent(in) :: bc_x, bc_y, bc_z
    real(8), dimension (:,:,:), intent(inout) :: u_star, v_star, w_star
    real(8), dimension (:,:), intent(in) :: bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz
    real(8), dimension (:,:), intent(in) :: bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz
    real(8), dimension (:,:), intent(in) :: bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz

    integer, dimension (3) :: nu, nv, nw, i

    nu=ubound(u_star)
    nv=ubound(v_star)
    nw=ubound(w_star)
    !!!!!!!!!!!!!! X component !!!!!!!!!!!!!!
    if (bc_x==1) then
        u_star(nu(1),:,:)=u_star(1      ,:,:)
        v_star(1    ,:,:)=v_star(nv(1)-1,:,:)
        v_star(nv(1),:,:)=v_star(2      ,:,:)
        w_star(1    ,:,:)=w_star(nw(1)-1,:,:)
        w_star(nw(1),:,:)=w_star(2      ,:,:)
    else if (bc_x==2) then
        u_star(1    ,:,:)=  bx_u_1
        u_star(nu(1),:,:)=  bx_u_nx
        v_star(1    ,:,:)=2*bx_v_1 -v_star(2      ,:,:)
        v_star(nv(1),:,:)=2*bx_v_nx-v_star(nv(1)-1,:,:)
        w_star(1    ,:,:)=2*bx_w_1 -w_star(2      ,:,:)
        w_star(nw(1),:,:)=2*bx_w_nx-w_star(nw(1)-1,:,:)
    else if (bc_x==3) then

    else if (bc_x==4) then
        u_star(1    ,:,:)=bx_u_1
        u_star(nu(1),:,:)=bx_u_nx
        v_star(1    ,:,:)=bx_v_1
        v_star(nv(1),:,:)=bx_v_nx
        w_star(1    ,:,:)=bx_w_1
        w_star(nw(1),:,:)=bx_w_nx
    end if

    !!!!!!!!!!!!!! Y component !!!!!!!!!!!!!!
    if (bc_y==1) then
        u_star(:,1    ,:)=u_star(:,nu(2)-1,:)
        u_star(:,nu(2),:)=u_star(:,2      ,:)
        v_star(:,nv(2),:)=v_star(:,1      ,:)
        w_star(:,1    ,:)=w_star(:,nw(2)-1,:)
        w_star(:,nw(2),:)=w_star(:,2      ,:)
    else if (bc_y==2) then
        u_star(:,1    ,:)=2*by_u_1 -u_star(:,2      ,:)
        u_star(:,nu(2),:)=2*by_u_ny-u_star(:,nu(2)-1,:)
        v_star(:,1    ,:)=  by_v_1
        v_star(:,nv(2),:)=  by_v_ny
        w_star(:,1    ,:)=2*by_w_1 -w_star(:,2      ,:)
        w_star(:,nw(2),:)=2*by_w_ny-w_star(:,nw(2)-1,:)
    else if (bc_y==3) then

    else if (bc_y==4) then
        u_star(:,1    ,:)=by_u_1
        u_star(:,nu(2),:)=by_u_ny
        v_star(:,1    ,:)=by_v_1
        v_star(:,nv(2),:)=by_v_ny
        w_star(:,1    ,:)=by_w_1
        w_star(:,nw(2),:)=by_w_ny
    end if

    !!!!!!!!!!!!!! Z component !!!!!!!!!!!!!!
    if (bc_z==1) then
        u_star(:,:,1    )=u_star(:,:,nu(3)-1)
        u_star(:,:,nu(3))=u_star(:,:,2      )
        v_star(:,:,1    )=v_star(:,:,nv(3)-1)
        v_star(:,:,nv(3))=v_star(:,:,2      )
        w_star(:,:,nw(3))=w_star(:,:,1      )
    else if (bc_z==2) then
        u_star(:,:,1    )=2*bz_u_1 -u_star(:,:,2      )
        u_star(:,:,nu(3))=2*bz_u_nz-u_star(:,:,nu(3)-1)
        v_star(:,:,1    )=2*bz_v_1 -v_star(:,:,2      )
        v_star(:,:,nv(3))=2*bz_v_nz-v_star(:,:,nv(3)-1)
        w_star(:,:,1    )=  bz_w_1
        w_star(:,:,nw(3))=  bz_w_nz
    else if (bc_z==3) then

    else if (bc_z==4) then
        u_star(:,:,1    )=bz_u_1
        u_star(:,:,nu(3))=bz_u_nz
        v_star(:,:,1    )=bz_v_1
        v_star(:,:,nv(3))=bz_v_nz
        w_star(:,:,1    )=bz_w_1
        w_star(:,:,nw(3))=bz_w_nz
    end if

    end subroutine vel_bc_staggered


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Poisson_LHS_staggered(nxp, nyp, nzp, dx2, dy2, dz2, pbc_x, pbc_y, pbc_z, dx, dy, dz) result(LHS_poisson)

    implicit none

    integer, intent(in) :: nxp, nyp, nzp, pbc_x, pbc_y, pbc_z
    real(8), intent(in) :: dx2, dy2, dz2, dx, dy, dz
    !real(8), dimension (:,:,:), intent(out) :: LHS_poisson
    class(csr), allocatable :: LHS_poisson
    class(coo), allocatable :: LHS_poisson_coo

    integer :: nx, ny, nz, i, j, k
    real(8) :: idx2, idy2, idz2
    real(8) :: Ix(nxp,nxp), Iy(nyp,nyp), Iz(nzp,nzp)
    real(8) :: Ax(nxp,nxp), Ay(nyp,nyp), Az(nzp,nzp)

    integer :: temp2
    integer, dimension (:,:), allocatable :: tempx, tempy, tempz

    integer, dimension(:), allocatable :: row1, row2, row3, row4, row5, row6, col1, col2, col3, col4, col5, col6, row_tol
    real(8), dimension((nxp*nyp*2+nxp*nzp*2+nyp*nzp*2)*2+1) :: vv
    integer, dimension((nxp*nyp*2+nxp*nzp*2+nyp*nzp*2)*2+1) :: ii, jj

    type(coo), allocatable :: co1, co2, co3, co4, co5, co6
    type(coo), allocatable :: t1, t2, t3, t4, t5, t6

    integer :: info=-111, job(8)=0

    !system_clock
    REAL(8) :: system_clock_rate
    INTEGER :: c1,c2,cr,cm

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    system_clock_rate = REAL(cr)

    Ix=0.0d0
    Iy=0.0d0
    Iz=0.0d0
    Ax=0.0d0
    Ay=0.0d0
    Az=0.0d0

    CALL SYSTEM_CLOCK(c1)
    print *, "Poisson_LHS_staggered start..."
    do i=1,nxp
        Ix(i,i)=1.0d0

        Ax(i,i)=-2.0d0
        if (i<nxp) then
            Ax(i+1,i)=1.0d0
            Ax(i,i+1)=1.0d0
        end if
    end do
    do i=1,nyp
        Iy(i,i)=1.0d0

        Ay(i,i)=-2.0d0
        if (i<nyp) then
            Ay(i+1,i)=1.0d0
            Ay(i,i+1)=1.0d0
        end if
    end do
    do i=1,nzp
        Iz(i,i)=1.0d0

        Az(i,i)=-2.0d0
        if (i<nzp) then
            Az(i+1,i)=1.0d0
            Az(i,i+1)=1.0d0
        end if
    end do

    Ix(1,:)=0.0d0; Ix(nxp,:)=0.0d0;
    Iy(1,:)=0.0d0; Iy(nyp,:)=0.0d0;
    Iz(1,:)=0.0d0; Iz(nzp,:)=0.0d0;
    Ax(1,:)=0.0d0; Ax(nxp,:)=0.0d0;
    Ay(1,:)=0.0d0; Ay(nyp,:)=0.0d0;
    Az(1,:)=0.0d0; Az(nzp,:)=0.0d0;
    !call print_mat(Ix)
    !call print_mat(Ay)
    
    Ax=Ax/dx2; Ay=Ay/dy2; Az=Az/dz2;

    CALL SYSTEM_CLOCK(c2)
    print '("   Get Ax/Ay/Az: " (f6.4) " second")', (c2-c1)/system_clock_rate

    !LHS_poisson = -kron(kron(Iz,Iy),Ax)-kron(kron(Iz,Ay),Ix)-kron(kron(Az,Iy),Ix);

    co1=dencoo(Ix)
    co2=dencoo(Iy)
    co3=dencoo(Iz)
    co4=dencoo(Ax)
    co5=dencoo(Ay)
    co6=dencoo(Az)
    t1=kron(co3,co2) !kron(Iz,Iy)
    t2=kron(co3,co5) !kron(Iz,Ay)
    t3=kron(co6,co2) !kron(Az,Iy)
    t4=kron(t1,co4); !s1=t4%to_csr() !kron(kron(Iz,Iy),Ax)
    t5=kron(t2,co1); !s2=t5%to_csr() !kron(kron(Iz,Ay),Ix)
    t6=kron(t3,co1); !s3=t6%to_csr() !kron(kron(Az,Iy),Ix)

    if (t4%shape(1)/=t5%shape(1) .or. t4%shape(1)/=t6%shape(1) .or. t4%shape(2)/=t5%shape(2) .or. t4%shape(2)/=t6%shape(2)) then
        print *, "COO sparse kron incorrect. (subroutine Poisson_LHS_staggered)"
        stop (1)
    end if

    LHS_poisson=add( add(t4%to_csr(), t5%to_csr()), t6%to_csr() )
    LHS_poisson_coo=from_csr(LHS_poisson)
    CALL SYSTEM_CLOCK(c2)
    print '("   Kron product: " (f6.4) " second")', (c2-c1)/system_clock_rate

    vv=0.0d0; ii=0; jj=0;

    temp2=1
    tempy=repmat_1d_2d( [(i, i=1, nyp)], 2, nzp )
    tempz=repmat_1d_2d( [(i, i=1, nzp)], 1, nyp )

    if (pbc_x==1) then

        row1=[1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        col1=[nxp-1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row1)
            where (ii==row1(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row1(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=row1
        vv(temp2:temp2+nyp*nzp-1)=1.0d0
        temp2=temp2+nyp*nzp
        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=col1
        vv(temp2:temp2+nyp*nzp-1)=-1.0d0
        temp2=temp2+nyp*nzp

        row2=[nxp+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        col2=[2+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row2)
            where (ii==row2(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row2(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=row2
        vv(temp2:temp2+nyp*nzp-1)=1.0d0
        temp2=temp2+nyp*nzp
        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=col2
        vv(temp2:temp2+nyp*nzp-1)=-1.0d0
        temp2=temp2+nyp*nzp

    else if (pbc_x==2) then
    else if (pbc_x==3) then
    else if (pbc_x==4) then

    end if

    tempx=repmat_1d_2d( [(i, i=1, nxp)], 2, nzp )
    tempz=repmat_1d_2d( [(i, i=1, nzp)], 1, nxp )
    if (pbc_y==1) then

        row3=[tempx+(1-1)*nxp+(tempz-1)*nxp*nyp]
        col3=[tempx+(nyp-1-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row3)
            where (ii==row3(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row3(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=row3
        vv(temp2:temp2+nxp*nzp-1)=1.0d0
        temp2=temp2+nxp*nzp
        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=col3
        vv(temp2:temp2+nxp*nzp-1)=-1.0d0
        temp2=temp2+nxp*nzp

        row4=[tempx+(nyp-1)*nxp+(tempz-1)*nxp*nyp]
        col4=[tempx+(2-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row4)
            where (ii==row4(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row4(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=row4
        vv(temp2:temp2+nxp*nzp-1)=1.0d0
        temp2=temp2+nxp*nzp
        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=col4
        vv(temp2:temp2+nxp*nzp-1)=-1.0d0
        temp2=temp2+nxp*nzp

    else if (pbc_y==2) then
    else if (pbc_y==3) then
    else if (pbc_y==4) then

    end if

    tempx=repmat_1d_2d( [(i, i=1, nxp)], 2, nzp )
    tempy=repmat_1d_2d( [(i, i=1, nyp)], 1, nyp )
    if (pbc_y==1) then

        row5=[tempx+(tempy-1)*nxp+(1-1)*nxp*nyp]
        col5=[tempx+(tempy-1)*nxp+(nzp-1-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row5)
            where (ii==row5(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row5(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=row5
        vv(temp2:temp2+nxp*nyp-1)=1.0d0
        temp2=temp2+nxp*nyp
        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=col5
        vv(temp2:temp2+nxp*nyp-1)=-1.0d0
        temp2=temp2+nxp*nyp

        row6=[tempx+(tempy-1)*nxp+(nzp-1)*nxp*nyp]
        col6=[tempx+(tempy-1)*nxp+(2-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row6)
            where (ii==row6(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row6(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=row6
        vv(temp2:temp2+nxp*nyp-1)=1.0d0
        temp2=temp2+nxp*nyp
        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=col6
        vv(temp2:temp2+nxp*nyp-1)=-1.0d0
        temp2=temp2+nxp*nyp

    else if (pbc_y==2) then
    else if (pbc_y==3) then
    else if (pbc_y==4) then

    end if

    !row_tol=unique_sort([row1, row2, row3, row4, row5, row6])
    !do i=1,size(row_tol)
    !    where (LHS_poisson_coo%row==row_tol(i)) LHS_poisson_coo%value=0.0
    !end do

    co1=coo_init_or_clean(LHS_poisson_coo%value, LHS_poisson_coo%row, LHS_poisson_coo%col, LHS_poisson_coo%shape)
    co2=coo_init_or_clean(vv, ii, jj, LHS_poisson_coo%shape)

    LHS_poisson=add( co1%to_csr(), co2%to_csr() )
    CALL SYSTEM_CLOCK(c1)
    print '("   LHS_poisson completed: " (f6.4) " second")', (c2-c1)/system_clock_rate
    print *, "**************************************"

    end function Poisson_LHS_staggered


    end module NS_functions

