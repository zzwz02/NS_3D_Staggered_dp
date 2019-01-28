    include 'mkl_dfti.f90'
    include 'mkl_trig_transforms.f90'

    module NS_functions

    use mkl_trig_transforms
    use FD_functions
    use coo_mod
    use csr_mod
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
    !$omp parallel do
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
    !$omp end parallel do

    end subroutine TGV

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cal_conv(u, v, w, bc_x, bc_y, bc_z, dx, dy, dz, conv_x, conv_y, conv_z)
    implicit none

    integer, intent(in) :: bc_x, bc_y, bc_z
    real(8), intent(in) :: dx, dy, dz
    real(8), dimension (:,:,:), intent(in) :: u, v, w
    real(8), dimension (:,:,:), intent(out), optional :: conv_x, conv_y, conv_z
    real(8), dimension (:,:,:), allocatable :: ux, uy, uz
    real(8), dimension (:,:,:), allocatable :: vx, vy, vz
    real(8), dimension (:,:,:), allocatable :: wx, wy, wz
    integer :: nx, ny, nz

    nx=size(u,1)-1; ny=size(v,2)-1; nz=size(w,3)-1
    if (.not. allocated(uy)) then
        allocate(uy(nx+1,ny+1,nz+2), uz(nx+1,ny+2,nz+1), vx(nx+1,ny+1,nz+2), vz(nx+2,ny+1,nz+1), wx(nx+1,ny+2,nz+1), wy(nx+2,ny+1,nz+1))
        if (bc_x==1) then
            allocate( ux(nx+1,ny+2,nz+2) )
        else
            allocate( ux(nx+0,ny+2,nz+2) )
        end if

        if (bc_y==1) then
            allocate( vy(nx+2,ny+1,nz+2) )
        else
            allocate( vy(nx+2,ny+0,nz+2) )
        end if

        if (bc_z==1) then
            allocate( wz(nx+2,ny+2,nz+1) )
        else
            allocate( wz(nx+2,ny+2,nz+0) )
        end if
    end if

    ux=avg(u,1,bc_x); uy=avg(u,2,0);    uz=avg(u,3,0);
    vx=avg(v,1,0);    vy=avg(v,2,bc_y); vz=avg(v,3,0);
    wx=avg(w,1,0);    wy=avg(w,2,0);    wz=avg(w,3,bc_z);

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

    end subroutine cal_conv

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cal_diff(u, v, w, bc_x, bc_y, bc_z, dx2, dy2, dz2, diff_x, diff_y, diff_z)
    implicit none

    integer, intent(in) :: bc_x, bc_y, bc_z
    real(8), intent(in) :: dx2, dy2, dz2
    real(8), dimension (:,:,:), intent(in) :: u, v, w
    real(8), dimension (:,:,:), intent(out), optional :: diff_x, diff_y, diff_z

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

    end subroutine cal_diff

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

    integer, dimension (3) :: nu, nv, nw

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
    subroutine vel_bc_staggered_CN2(u_star, v_star, w_star, u_1, u_end, v_1, v_end, w_1, w_end, bc, dir)

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall)
    ! pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell

    implicit none
    integer, intent(in) :: bc, dir
    real(8), dimension (:,:,:), intent(inout) :: u_star, v_star, w_star
    real(8), dimension (:,:), intent(in) :: u_1, u_end, v_1, v_end, w_1, w_end

    integer, dimension (3) :: nu, nv, nw

    nu=ubound(u_star)
    nv=ubound(v_star)
    nw=ubound(w_star)

    !!!!!!!!!!!!!! X component !!!!!!!!!!!!!!
    if (dir==1) then
        if (bc==1) then
            u_star(nu(1),:,:)=0.0d0
            v_star(1    ,:,:)=0.0d0
            v_star(nv(1),:,:)=0.0d0
            w_star(1    ,:,:)=0.0d0
            w_star(nw(1),:,:)=0.0d0
        else if (bc==2) then
            u_star(1    ,:,:)=u_1
            u_star(nu(1),:,:)=u_end
            v_star(1    ,:,:)=v_1
            v_star(nv(1),:,:)=v_end
            w_star(1    ,:,:)=w_1
            w_star(nw(1),:,:)=w_end
        else if (bc==3) then

        else if (bc==4) then
            u_star(1    ,:,:)=u_1
            u_star(nu(1),:,:)=u_end
            v_star(1    ,:,:)=v_1
            v_star(nv(1),:,:)=v_end
            w_star(1    ,:,:)=w_1
            w_star(nw(1),:,:)=w_end
        end if

        !!!!!!!!!!!!!! Y component !!!!!!!!!!!!!!
    else if (dir==2) then
        if (bc==1) then
            u_star(:,1    ,:)=0.0d0
            u_star(:,nu(2),:)=0.0d0
            v_star(:,nv(2),:)=0.0d0
            w_star(:,1    ,:)=0.0d0
            w_star(:,nw(2),:)=0.0d0
        else if (bc==2) then
            u_star(:,1    ,:)=u_1
            u_star(:,nu(2),:)=u_end
            v_star(:,1    ,:)=v_1
            v_star(:,nv(2),:)=v_end
            w_star(:,1    ,:)=w_1
            w_star(:,nw(2),:)=w_end
        else if (bc==3) then

        else if (bc==4) then
            u_star(:,1    ,:)=u_1
            u_star(:,nu(2),:)=u_end
            v_star(:,1    ,:)=v_1
            v_star(:,nv(2),:)=v_end
            w_star(:,1    ,:)=w_1
            w_star(:,nw(2),:)=w_end
        end if

        !!!!!!!!!!!!!! Z component !!!!!!!!!!!!!!
    else if (dir==3) then
        if (bc==1) then
            u_star(:,:,1    )=0.0d0
            u_star(:,:,nu(3))=0.0d0
            v_star(:,:,1    )=0.0d0
            v_star(:,:,nv(3))=0.0d0
            w_star(:,:,nw(3))=0.0d0
        else if (bc==2) then
            u_star(:,:,1    )=u_1
            u_star(:,:,nu(3))=u_end
            v_star(:,:,1    )=v_1
            v_star(:,:,nv(3))=v_end
            w_star(:,:,1    )=w_1
            w_star(:,:,nw(3))=w_end
        else if (bc==3) then

        else if (bc==4) then
            u_star(:,:,1    )=u_1
            u_star(:,:,nu(3))=u_end
            v_star(:,:,1    )=v_1
            v_star(:,:,nv(3))=v_end
            w_star(:,:,1    )=w_1
            w_star(:,:,nw(3))=w_end
        end if
    else
        print *, "dir must be 1/2/3. (SUBROUTINE vel_bc_staggered_CN2)"
    end if

    end subroutine vel_bc_staggered_CN2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_velpr_bc(u_bc, v_bc, w_bc, p_bc, bc_x, bc_y, bc_z, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, &
        bx_u_1, bx_u_nx, by_u_1, by_u_ny, bz_u_1, bz_u_nz, bx_v_1, bx_v_nx, by_v_1, by_v_ny, bz_v_1, bz_v_nz, &
        bx_w_1, bx_w_nx, by_w_1, by_w_ny, bz_w_1, bz_w_nz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)
    implicit none

    integer, intent(in) :: bc_x, bc_y, bc_z, pbc_x, pbc_y, pbc_z, nx, ny, nz
    real(8), intent(in) :: dx, dy, dz
    real(8), intent(in) :: u_bc(nx+1, ny+2, nz+2), v_bc(nx+2, ny+1, nz+2), w_bc(nx+2, ny+2, nz+1), p_bc(nx+2, ny+2, nz+2)
    real(8), intent(out) :: bx_u_1(ny+2,nz+2), bx_u_nx(ny+2,nz+2), by_u_1(nx+1,nz+2), by_u_ny(nx+1,nz+2), bz_u_1(nx+1,ny+2), bz_u_nz(nx+1,ny+2)
    real(8), intent(out) :: bx_v_1(ny+1,nz+2), bx_v_nx(ny+1,nz+2), by_v_1(nx+2,nz+2), by_v_ny(nx+2,nz+2), bz_v_1(nx+2,ny+1), bz_v_nz(nx+2,ny+1)
    real(8), intent(out) :: bx_w_1(ny+2,nz+1), bx_w_nx(ny+2,nz+1), by_w_1(nx+2,nz+1), by_w_ny(nx+2,nz+1), bz_w_1(nx+2,ny+2), bz_w_nz(nx+2,ny+2)
    real(8), intent(out) :: bx_p_1(ny+2,nz+2), bx_p_nx(ny+2,nz+2), by_p_1(nx+2,nz+2), by_p_ny(nx+2,nz+2), bz_p_1(nx+2,ny+2), bz_p_nz(nx+2,ny+2)

    integer :: nxp, nyp, nzp

    nxp=nx+2; nyp=ny+2; nzp=nz+2

    if (bc_x==2) then
        bx_u_1 =u_bc(1,:,:);
        bx_u_nx=u_bc(nx+1,:,:);
        by_u_1 =(u_bc(:,1,:)+u_bc(:,2,:))/2;
        by_u_ny=(u_bc(:,ny+2,:)+u_bc(:,ny+2-1,:))/2;
        bz_u_1 =(u_bc(:,:,1)+u_bc(:,:,2))/2;
        bz_u_nz=(u_bc(:,:,nz+2)+u_bc(:,:,nz+2-1))/2;

        bx_v_1 =(v_bc(1,:,:)+v_bc(2,:,:))/2;
        bx_v_nx=(v_bc(nx+2,:,:)+v_bc(nx+2-1,:,:))/2;
        by_v_1 =v_bc(:,1,:);
        by_v_ny=v_bc(:,ny+1,:);
        bz_v_1 =(v_bc(:,:,1)+v_bc(:,:,2))/2;
        bz_v_nz=(v_bc(:,:,nz+2)+v_bc(:,:,nz+2-1))/2;

        bx_w_1 =(w_bc(1,:,:)+w_bc(2,:,:))/2;
        bx_w_nx=(w_bc(nx+2,:,:)+w_bc(nx+2-1,:,:))/2;
        by_w_1 =(w_bc(:,1,:)+w_bc(:,2,:))/2;
        by_w_ny=(w_bc(:,ny+2,:)+w_bc(:,ny+2-1,:))/2;
        bz_w_1 =w_bc(:,:,1);
        bz_w_nz=w_bc(:,:,nz+1);
    else if (bc_x==4) then
        bx_u_1 =u_bc(1,:,:);
        bx_u_nx=u_bc(nx+1,:,:);
        by_u_1 =u_bc(:,1,:);
        by_u_ny=u_bc(:,ny+2,:);
        bz_u_1 =u_bc(:,:,1);
        bz_u_nz=u_bc(:,:,nz+2);

        bx_v_1 =v_bc(1,:,:);
        bx_v_nx=v_bc(nx+2,:,:);
        by_v_1 =v_bc(:,1,:);
        by_v_ny=v_bc(:,ny+1,:);
        bz_v_1 =v_bc(:,:,1);
        bz_v_nz=v_bc(:,:,nz+2);

        bx_w_1 =w_bc(1,:,:);
        bx_w_nx=w_bc(nx+2,:,:);
        by_w_1 =w_bc(:,1,:);
        by_w_ny=w_bc(:,ny+2,:);
        bz_w_1 =w_bc(:,:,1);
        bz_w_nz=w_bc(:,:,nz+1);
    end if

    call get_pr_bc(p_bc, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)

    end subroutine get_velpr_bc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine get_pr_bc(p_bc, pbc_x, pbc_y, pbc_z, nx, ny, nz, dx, dy, dz, bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz)
    implicit none

    integer, intent(in) :: pbc_x, pbc_y, pbc_z, nx, ny, nz
    real(8), intent(in) :: dx, dy, dz
    real(8), dimension (nx+2,ny+2,nz+2), intent(in) :: p_bc
    real(8), intent(out) :: bx_p_1(ny+2,nz+2), bx_p_nx(ny+2,nz+2), by_p_1(nx+2,nz+2), by_p_ny(nx+2,nz+2), bz_p_1(nx+2,ny+2), bz_p_nz(nx+2,ny+2)

    integer :: nxp, nyp, nzp

    nxp=nx+2; nyp=ny+2; nzp=nz+2

    if (pbc_x==2) then
        bx_p_1 =(p_bc(1,:,:)  +p_bc(2,:,:))/2;
        bx_p_nx=(p_bc(nxp,:,:)+p_bc(nxp-1,:,:))/2;
    else if (pbc_x==4) then
        bx_p_1 =p_bc(1,:,:);
        bx_p_nx=p_bc(nxp,:,:);
    else if (pbc_x==3) then
        bx_p_1 =reshape([diff(p_bc(1:2,:,:),1,1)/dx]      , [nyp,nzp])
        bx_p_nx=reshape([diff(p_bc(nxp-1:nxp,:,:),1,1)/dx], [nyp,nzp])

    end if
    
    if (pbc_y==2) then
        by_p_1 =(p_bc(:,1,:)  +p_bc(:,2,:))/2;
        by_p_ny=(p_bc(:,nyp,:)+p_bc(:,nyp-1,:))/2;
    else if (pbc_y==4) then
        by_p_1 =p_bc(:,1,:);
        by_p_ny=p_bc(:,nyp,:);
    else if (pbc_y==3) then
        by_p_1 =reshape([diff(p_bc(:,1:2,:),1,3)/dy]      , [nxp,nzp])
        by_p_ny=reshape([diff(p_bc(:,nyp-1:nyp,:),1,3)/dy], [nxp,nzp])
    end if
    
    if (pbc_z==2) then
        bz_p_1 =(p_bc(:,:,1)  +p_bc(:,:,2))/2;
        bz_p_nz=(p_bc(:,:,nzp)+p_bc(:,:,nzp-1))/2;
    else if (pbc_z==4) then
        bz_p_1 =p_bc(:,:,1);
        bz_p_nz=p_bc(:,:,nzp);
    else if (pbc_z==3) then
        bz_p_1 =reshape([diff(p_bc(:,:,1:2),1,3)/dz]      , [nxp,nyp])
        bz_p_nz=reshape([diff(p_bc(:,:,nzp-1:nzp),1,3)/dz], [nxp,nyp])
    end if

    end subroutine get_pr_bc

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
    print '("   Get Ax/Ay/Az: ", F8.4, " second")', (c2-c1)/system_clock_rate

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
    print '("   Kron product: ", F8.4, " second")', (c2-c1)/system_clock_rate

    vv=0.0d0; ii=0; jj=0;

    temp2=1
    tempy=spread( [(i, i=1, nyp)], 2, nzp )
    tempz=spread( [(i, i=1, nzp)], 1, nyp )
    if (pbc_x==1) then ! Periodic

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

    else if (pbc_x==2) then ! Dirichlet on boundary (cell wall)

        row1=[1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        col1=[2+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row1)
            where (ii==row1(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row1(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=row1
        vv(temp2:temp2+nyp*nzp-1)=0.5d0
        temp2=temp2+nyp*nzp
        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=col1
        vv(temp2:temp2+nyp*nzp-1)=0.5d0
        temp2=temp2+nyp*nzp

        row2=[nxp+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        col2=[nxp-1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row2)
            where (ii==row2(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row2(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=row2
        vv(temp2:temp2+nyp*nzp-1)=0.5d0
        temp2=temp2+nyp*nzp
        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=col2
        vv(temp2:temp2+nyp*nzp-1)=0.5d0
        temp2=temp2+nyp*nzp

    else if (pbc_x==3) then ! Neumann on boundary (cell wall)

        row1=[1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        col1=[2+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row1)
            where (ii==row1(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row1(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=row1
        vv(temp2:temp2+nyp*nzp-1)=-1.0d0/dx
        temp2=temp2+nyp*nzp
        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=col1
        vv(temp2:temp2+nyp*nzp-1)=1.0d0/dx
        temp2=temp2+nyp*nzp

        row2=[nxp+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        col2=[nxp-1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row2)
            where (ii==row2(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row2(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=row2
        vv(temp2:temp2+nyp*nzp-1)=1.0d0/dx
        temp2=temp2+nyp*nzp
        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=col2
        vv(temp2:temp2+nyp*nzp-1)=-1.0d0/dx
        temp2=temp2+nyp*nzp

    else if (pbc_x==4) then ! Dirichlet on ghost cell

        row1=[1+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row1)
            where (ii==row1(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row1(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row1
        jj(temp2:temp2+nyp*nzp-1)=row1
        vv(temp2:temp2+nyp*nzp-1)=1.0d0
        temp2=temp2+nyp*nzp

        row2=[nxp+(tempy-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row2)
            where (ii==row2(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row2(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nyp*nzp-1)=row2
        jj(temp2:temp2+nyp*nzp-1)=row2
        vv(temp2:temp2+nyp*nzp-1)=1.0d0
        temp2=temp2+nyp*nzp

    end if

    tempx=spread( [(i, i=1, nxp)], 2, nzp )
    tempz=spread( [(i, i=1, nzp)], 1, nxp )
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

        row3=[tempx+(1-1)*nxp+(tempz-1)*nxp*nyp]
        col3=[tempx+(2-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row3)
            where (ii==row3(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row3(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=row3
        vv(temp2:temp2+nxp*nzp-1)=0.5d0
        temp2=temp2+nxp*nzp
        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=col3
        vv(temp2:temp2+nxp*nzp-1)=0.5d0
        temp2=temp2+nxp*nzp

        row4=[tempx+(nyp-1)*nxp+(tempz-1)*nxp*nyp]
        col4=[tempx+(nyp-1-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row4)
            where (ii==row4(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row4(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=row4
        vv(temp2:temp2+nxp*nzp-1)=0.5d0
        temp2=temp2+nxp*nzp
        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=col4
        vv(temp2:temp2+nxp*nzp-1)=0.5d0
        temp2=temp2+nxp*nzp

    else if (pbc_y==3) then

        row3=[tempx+(1-1)*nxp+(tempz-1)*nxp*nyp]
        col3=[tempx+(2-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row3)
            where (ii==row3(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row3(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=row3
        vv(temp2:temp2+nxp*nzp-1)=-1.0d0/dy
        temp2=temp2+nxp*nzp
        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=col3
        vv(temp2:temp2+nxp*nzp-1)=1.0d0/dy
        temp2=temp2+nxp*nzp

        row4=[tempx+(nyp-1)*nxp+(tempz-1)*nxp*nyp]
        col4=[tempx+(nyp-1-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row4)
            where (ii==row4(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row4(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=row4
        vv(temp2:temp2+nxp*nzp-1)=1.0d0/dy
        temp2=temp2+nxp*nzp
        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=col4
        vv(temp2:temp2+nxp*nzp-1)=-1.0d0/dy
        temp2=temp2+nxp*nzp

    else if (pbc_y==4) then

        row3=[tempx+(1-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row3)
            where (ii==row3(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row3(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row3
        jj(temp2:temp2+nxp*nzp-1)=row3
        vv(temp2:temp2+nxp*nzp-1)=1.0d0
        temp2=temp2+nxp*nzp

        row4=[tempx+(nyp-1)*nxp+(tempz-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row4)
            where (ii==row4(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row4(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nzp-1)=row4
        jj(temp2:temp2+nxp*nzp-1)=row4
        vv(temp2:temp2+nxp*nzp-1)=1.0d0
        temp2=temp2+nxp*nzp

    end if

    tempx=spread( [(i, i=1, nxp)], 2, nzp )
    tempy=spread( [(i, i=1, nyp)], 1, nyp )
    if (pbc_z==1) then

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

    else if (pbc_z==2) then

        row5=[tempx+(tempy-1)*nxp+(1-1)*nxp*nyp]
        col5=[tempx+(tempy-1)*nxp+(2-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row5)
            where (ii==row5(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row5(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=row5
        vv(temp2:temp2+nxp*nyp-1)=0.5d0
        temp2=temp2+nxp*nyp
        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=col5
        vv(temp2:temp2+nxp*nyp-1)=0.5d0
        temp2=temp2+nxp*nyp

        row6=[tempx+(tempy-1)*nxp+(nzp-1)*nxp*nyp]
        col6=[tempx+(tempy-1)*nxp+(nzp-1-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row6)
            where (ii==row6(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row6(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=row6
        vv(temp2:temp2+nxp*nyp-1)=0.5d0
        temp2=temp2+nxp*nyp
        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=col6
        vv(temp2:temp2+nxp*nyp-1)=0.5d0
        temp2=temp2+nxp*nyp

    else if (pbc_z==3) then

        row5=[tempx+(tempy-1)*nxp+(1-1)*nxp*nyp]
        col5=[tempx+(tempy-1)*nxp+(2-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row5)
            where (ii==row5(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row5(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=row5
        vv(temp2:temp2+nxp*nyp-1)=-1.0d0/dz
        temp2=temp2+nxp*nyp
        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=col5
        vv(temp2:temp2+nxp*nyp-1)=1.0d0/dz
        temp2=temp2+nxp*nyp

        row6=[tempx+(tempy-1)*nxp+(nzp-1)*nxp*nyp]
        col6=[tempx+(tempy-1)*nxp+(nzp-1-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row6)
            where (ii==row6(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row6(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=row6
        vv(temp2:temp2+nxp*nyp-1)=1.0d0/dz
        temp2=temp2+nxp*nyp
        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=col6
        vv(temp2:temp2+nxp*nyp-1)=-1.0d0/dz
        temp2=temp2+nxp*nyp

    else if (pbc_z==4) then

        row5=[tempx+(tempy-1)*nxp+(1-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row5)
            where (ii==row5(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row5(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row5
        jj(temp2:temp2+nxp*nyp-1)=row5
        vv(temp2:temp2+nxp*nyp-1)=1.0d0
        temp2=temp2+nxp*nyp

        row6=[tempx+(tempy-1)*nxp+(nzp-1)*nxp*nyp]
        !call print_mat(float(row))
        do i=1,size(row6)
            where (ii==row6(i)) vv=0.0d0
            !where (LHS_poisson_coo%row==row6(i)) LHS_poisson_coo%value=0.0
        end do

        ii(temp2:temp2+nxp*nyp-1)=row6
        jj(temp2:temp2+nxp*nyp-1)=row6
        vv(temp2:temp2+nxp*nyp-1)=1.0d0
        temp2=temp2+nxp*nyp

    end if

    !row_tol=unique_sort([row1, row2, row3, row4, row5, row6])
    !do i=1,size(row_tol)
    !    where (LHS_poisson_coo%row==row_tol(i)) LHS_poisson_coo%value=0.0
    !end do

    co1=coo_init_or_clean(LHS_poisson_coo%value, LHS_poisson_coo%row, LHS_poisson_coo%col, LHS_poisson_coo%shape)
    co2=coo_init_or_clean(vv, ii, jj, LHS_poisson_coo%shape)

    LHS_poisson=add( co1%to_csr(), co2%to_csr() )
    CALL SYSTEM_CLOCK(c2)
    print '("   LHS_poisson completed: ", F8.4, " second")', (c2-c1)/system_clock_rate
    print *, "**************************************"

    end function Poisson_LHS_staggered

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mat_CN(nx, bc_x, dx2, dx, A_mat, B_mat, dt0, nu)

    implicit none

    integer, intent(in) :: nx, bc_x
    real(8), intent(in) :: dx2, dx, dt0, nu
    real(8), intent(out) :: A_mat(nx+1,nx+1), B_mat(nx+2,nx+2)

    integer :: i
    real(8) :: idx2, idx

    idx2=1.0d0/dx2; idx=1.0/dx;

    A_mat=0; B_mat=0;

    do i=2,nx
        A_mat(i,i)=1.0d0+0.5d0*dt0*nu*2.0d0*idx2
        A_mat(i,i-1)=-0.5d0*dt0*nu*1.0d0*idx2
        A_mat(i,i+1)=-0.5d0*dt0*nu*1.0d0*idx2
    end do

    do i=2,nx+1
        B_mat(i,i)=1.0d0+0.5d0*dt0*nu*2.0d0*idx2;
        B_mat(i,i-1)=-0.5d0*dt0*nu*1.0d0*idx2;
        B_mat(i,i+1)=-0.5d0*dt0*nu*1.0d0*idx2;
    end do

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall); pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    if (bc_x==1) then
        A_mat(1,1)=1.0d0+0.5d0*dt0*nu*2.0d0*idx2; A_mat(1,2)=-0.5d0*dt0*nu*1.0d0*idx2; A_mat(1,nx+1)=-0.5d0*dt0*nu*1.0d0*idx2;
        A_mat(nx+1,nx+1)=1.0d0; A_mat(nx+1,1)=-1.0d0;

        B_mat(1,1)=1.0d0;       B_mat(1,nx+1)=-1.0d0;
        B_mat(nx+2,nx+2)=1.0d0; B_mat(nx+2,2)=-1.0d0;

    else if (bc_x==2) then
        A_mat(1,1)=1.0d0;       A_mat(nx+1,nx+1)=1.0d0;

        B_mat(1,1)=0.5d0;       B_mat(1,2)=0.5d0;
        B_mat(nx+2,nx+2)=0.5d0; B_mat(nx+2,nx+1)=0.5d0;

    else if (bc_x==3) then

    else if (bc_x==4) then
        A_mat(1,1)=1.0d0;       A_mat(nx+1,nx+1)=1.0d0;

        B_mat(1,1)=1.0d0;       B_mat(nx+2,nx+2)=1.0d0;
    end if

    end subroutine mat_CN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mat_CN_tridiagonal(nx, bc_x, dx2, dx, A_low, A_d, A_up, B_low, B_d, B_up, dt0, nu)

    implicit none

    integer, intent(in) :: nx, bc_x
    real(8), intent(in) :: dx2, dx, dt0, nu
    real(8), intent(out) :: A_low(nx), A_d(nx+1), A_up(nx), B_low(nx+1), B_d(nx+2), B_up(nx+1)

    integer :: i
    real(8) :: idx2, idx

    idx2=1.0d0/dx2; idx=1.0/dx;

    do i=2,nx
        A_d(i)=1.0d0+0.5d0*dt0*nu*2.0d0*idx2
        A_low(i-1)= -0.5d0*dt0*nu*1.0d0*idx2
        A_up(i)=    -0.5d0*dt0*nu*1.0d0*idx2
    end do

    do i=2,nx+1
        B_d(i)=1.0d0+0.5d0*dt0*nu*2.0d0*idx2;
        B_low(i-1)= -0.5d0*dt0*nu*1.0d0*idx2;
        B_up(i)=    -0.5d0*dt0*nu*1.0d0*idx2;
    end do

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall); pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell
    if (bc_x==1) then


    else if (bc_x==2) then
        A_d(1)=1.0d0;    A_up(1)=0.0d0;
        A_d(nx+1)=1.0d0; A_low(nx)=0.0d0;

        B_d(1)=0.5d0;    B_up(1)=0.5d0;
        B_d(nx+2)=0.5d0; B_low(nx+1)=0.5d0;
    else if (bc_x==3) then

    else if (bc_x==4) then
        A_d(1)=1.0d0;    A_up(1)=0.0d0;
        A_d(nx+1)=1.0d0; A_low(nx)=0.0d0;

        B_d(1)=1.0d0;    B_up(1)=0.0d0;
        B_d(nx+2)=1.0d0; B_low(nx+1)=0.0d0;
    end if

    end subroutine mat_CN_tridiagonal

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pr_bc_staggered_modify_rhs(RHS_poisson_internal, &
        bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz,&
        pbc_x, pbc_y, pbc_z, dx, dy, dz, dx2, dy2, dz2)

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall)
    ! pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell

    implicit none
    integer, intent(in) :: pbc_x, pbc_y, pbc_z
    real(8), dimension (:,:,:), intent(inout) :: RHS_poisson_internal
    real(8), dimension (:,:), intent(in) :: bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz
    real(8), intent(in) :: dx, dy, dz, dx2, dy2, dz2

    integer, dimension (3) :: np

    np=ubound(RHS_poisson_internal)
    !!!!!!!!!!!!!! X component !!!!!!!!!!!!!!
    if (pbc_x==1) then
    else if (pbc_x==2) then
        RHS_poisson_internal(1    ,:,:)=RHS_poisson_internal(1    ,:,:)-2*bx_p_1/dx2
        RHS_poisson_internal(np(1),:,:)=RHS_poisson_internal(np(1),:,:)-2*bx_p_nx/dx2
    else if (pbc_x==3) then
        RHS_poisson_internal(1    ,:,:)=RHS_poisson_internal(1    ,:,:)+bx_p_1/dx
        RHS_poisson_internal(np(1),:,:)=RHS_poisson_internal(np(1),:,:)-bx_p_nx/dx
    else if (pbc_x==4) then
        RHS_poisson_internal(1    ,:,:)=RHS_poisson_internal(1    ,:,:)-bx_p_1/dx2
        RHS_poisson_internal(np(1),:,:)=RHS_poisson_internal(np(1),:,:)-bx_p_nx/dx2
    end if

    !!!!!!!!!!!!!! Y component !!!!!!!!!!!!!!
    if (pbc_y==1) then
    else if (pbc_y==2) then
        RHS_poisson_internal(:,1    ,:)=RHS_poisson_internal(:,1    ,:)-2*by_p_1/dy2
        RHS_poisson_internal(:,np(2),:)=RHS_poisson_internal(:,np(2),:)-2*by_p_ny/dy2
    else if (pbc_y==3) then
        RHS_poisson_internal(:,1    ,:)=RHS_poisson_internal(:,1    ,:)+by_p_1/dy
        RHS_poisson_internal(:,np(2),:)=RHS_poisson_internal(:,np(2),:)-by_p_ny/dy
    else if (pbc_y==4) then
        RHS_poisson_internal(:,1    ,:)=RHS_poisson_internal(:,1    ,:)-by_p_1/dy2
        RHS_poisson_internal(:,np(2),:)=RHS_poisson_internal(:,np(2),:)-by_p_ny/dy2
    end if

    !!!!!!!!!!!!!! Z component !!!!!!!!!!!!!!
    if (pbc_z==1) then
    else if (pbc_z==2) then
        RHS_poisson_internal(:,:,1    )=RHS_poisson_internal(:,:,1    )-2*bz_p_1/dz2
        RHS_poisson_internal(:,:,np(3))=RHS_poisson_internal(:,:,np(3))-2*bz_p_nz/dz2
    else if (pbc_z==3) then
        RHS_poisson_internal(:,:,1    )=RHS_poisson_internal(:,:,1    )+bz_p_1/dz
        RHS_poisson_internal(:,:,np(3))=RHS_poisson_internal(:,:,np(3))-bz_p_nz/dz
    else if (pbc_z==4) then
        RHS_poisson_internal(:,:,1    )=RHS_poisson_internal(:,:,1    )-bz_p_1/dz2
        RHS_poisson_internal(:,:,np(3))=RHS_poisson_internal(:,:,np(3))-bz_p_nz/dz2
    end if

    end subroutine pr_bc_staggered_modify_rhs

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pr_bc_staggered(dp, &
        bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz,&
        pbc_x, pbc_y, pbc_z, dx, dy, dz)

    ! pbc=1 Periodic; pbc=2 Dirichlet on boundary (cell wall)
    ! pbc=3 Neumann on boundary (cell wall); pbc=4 Dirichlet on ghost cell

    implicit none
    integer, intent(in) :: pbc_x, pbc_y, pbc_z
    real(8), dimension (:,:,:), intent(inout) :: dp
    real(8), dimension (:,:), intent(in) :: bx_p_1, bx_p_nx, by_p_1, by_p_ny, bz_p_1, bz_p_nz
    real(8), intent(in) :: dx, dy, dz

    integer, dimension (3) :: np

    np=ubound(dp)
    !!!!!!!!!!!!!! X component !!!!!!!!!!!!!!
    if (pbc_x==1) then
        dp(1    ,:,:)=dp(np(1)-1,:,:)
        dp(np(1),:,:)=dp(2      ,:,:)
    else if (pbc_x==2) then
        dp(1    ,:,:)=2*bx_p_1 -dp(2      ,:,:)
        dp(np(1),:,:)=2*bx_p_nx-dp(np(1)-1,:,:)
    else if (pbc_x==3) then
        dp(1    ,:,:)=-dx*bx_p_1 +dp(2      ,:,:)
        dp(np(1),:,:)= dx*bx_p_nx+dp(np(1)-1,:,:)
    else if (pbc_x==4) then
        dp(1    ,:,:)=bx_p_1
        dp(np(1),:,:)=bx_p_nx
    end if

    !!!!!!!!!!!!!! Y component !!!!!!!!!!!!!!
    if (pbc_y==1) then
        dp(:,1    ,:)=dp(:,np(2)-1,:)
        dp(:,np(2),:)=dp(:,2      ,:)
    else if (pbc_y==2) then
        dp(:,1    ,:)=2*by_p_1 -dp(:,2      ,:)
        dp(:,np(2),:)=2*by_p_ny-dp(:,np(2)-1,:)
    else if (pbc_y==3) then
        dp(:,1    ,:)=-dy*by_p_1 +dp(:,2      ,:)
        dp(:,np(2),:)= dy*by_p_ny+dp(:,np(2)-1,:)
    else if (pbc_y==4) then
        dp(:,1    ,:)=by_p_1
        dp(:,np(2),:)=by_p_ny
    end if

    !!!!!!!!!!!!!! Z component !!!!!!!!!!!!!!
    if (pbc_z==1) then
        dp(:,:,1    )=dp(:,:,np(3)-1)
        dp(:,:,np(3))=dp(:,:,2      )
    else if (pbc_z==2) then
        dp(:,:,1    )=2*bz_p_1 -dp(:,:,2      )
        dp(:,:,np(3))=2*bz_p_nz-dp(:,:,np(3)-1)
    else if (pbc_z==3) then
        dp(:,:,1    )=-dz*bz_p_1 +dp(:,:,2      )
        dp(:,:,np(3))= dz*bz_p_nz+dp(:,:,np(3)-1)
    else if (pbc_z==4) then
        dp(:,:,1    )=bz_p_1
        dp(:,:,np(3))=bz_p_nz
    end if

    end subroutine pr_bc_staggered

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine DST_DCT_2(f, nx, handle, ipar, dpar, DST2, DCT2, forward, backward)

    implicit none

    logical, intent(in), optional :: DST2, DCT2, forward, backward
    real(8), intent(inout) :: f(nx+1), dpar(5*nx/2+2)
    !!!!!!!!!!!!!!! INTEL mkl_tt !!!!!!!!!!!!!!!
    integer, intent(inout) :: ipar(128)!, direction
    type(dfti_descriptor), pointer, intent(in) :: handle

    integer :: nx, status

    if (present(DST2) .and. DST2) then

        if (present(forward) .and. forward) then
            f(1:nx:2)=-f(1:nx:2)
            call d_backward_trig_transform(f,handle,ipar,dpar,status)
            f(1:nx)=f(nx:1:-1)
        elseif (present(backward) .and. backward) then
            f(1:nx)=f(nx:1:-1)
            call d_forward_trig_transform(f,handle,ipar,dpar,status)
            f(1:nx:2)=-f(1:nx:2)
        end if

    elseif (present(DCT2).and. DCT2) then

        if (present(forward) .and. forward) then
            call d_backward_trig_transform(f,handle,ipar,dpar,status)
        elseif (present(backward) .and. backward) then
            call d_forward_trig_transform(f,handle,ipar,dpar,status)
        end if

    end if

    end subroutine DST_DCT_2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine DCT_poisson_solver(rhs_tt, eig_tt, handle_x, handle_y, handle_z, &
        ipar_x, ipar_y, ipar_z, dpar_x, dpar_y, dpar_z, nx, ny, nz, pbc_x, pbc_y, pbc_z)

    implicit none

    integer, intent(in) :: nx, ny, nz, pbc_x, pbc_y, pbc_z
    real(8), intent(in) :: eig_tt(nx,ny,nz), dpar_x(5*nx/2+2), dpar_y(5*nx/2+2), dpar_z(5*nz/2+2)
    real(8), intent(out) :: rhs_tt(nx+1,ny+1,nz+1)
    integer, intent(inout) :: ipar_x(128), ipar_y(128), ipar_z(128)
    type(dfti_descriptor), pointer, intent(in) :: handle_x, handle_y, handle_z

    integer :: status, i, j, k

    !x-direction
    do k=1,nz
        do j=1,ny
            if (pbc_x==2) then ! DST-2 = flip sign of every other input, iDCT-3, reverse order of output
                rhs_tt(1:nx:2,j,k)=-rhs_tt(1:nx:2,j,k)
                call d_backward_trig_transform(rhs_tt(:,j,k),handle_x,ipar_x,dpar_x,status)
                rhs_tt(1:nx,j,k)=rhs_tt(nx:1:-1,j,k)
            elseif (pbc_x==3) then ! DCT-2 = iDCT-3
                call d_backward_trig_transform(rhs_tt(:,j,k),handle_x,ipar_x,dpar_x,status)
            elseif (pbc_x==4) then ! DST-1
                call d_forward_trig_transform(rhs_tt(:,j,k),handle_x,ipar_x,dpar_x,status)
            end if
        end do
    end do

    !y-direction
    do k=1,nz
        do i=1,nx
            if (pbc_y==2) then ! DST-2 = flip sign of every other input, iDCT-3, reverse order of output
                rhs_tt(i,1:ny:2,k)=-rhs_tt(i,1:ny:2,k)
                call d_backward_trig_transform(rhs_tt(i,:,k),handle_y,ipar_y,dpar_y,status)
                rhs_tt(i,1:ny,k)=rhs_tt(i,ny:1:-1,k)
            elseif (pbc_y==3) then ! DCT-2 = iDCT-3
                call d_backward_trig_transform(rhs_tt(i,:,k),handle_y,ipar_y,dpar_y,status)
            elseif (pbc_y==4) then ! DST-1
                call d_forward_trig_transform(rhs_tt(i,:,k),handle_y,ipar_y,dpar_y,status)
            end if
        end do
    end do

    !z-direction
    do j=1,ny
        do i=1,nx
            if (pbc_z==2) then ! DST-2 = flip sign of every other input, iDCT-3, reverse order of output
                rhs_tt(i,j,1:nz:2)=-rhs_tt(i,j,1:nz:2)
                call d_backward_trig_transform(rhs_tt(i,j,:),handle_z,ipar_z,dpar_z,status)
                !rhs_tt(i,j,1:nz)=rhs_tt(i,j,nz:1:-1)
            elseif (pbc_z==3) then ! DCT-2 = iDCT-3
                call d_backward_trig_transform(rhs_tt(i,j,:),handle_z,ipar_z,dpar_z,status)
            elseif (pbc_z==4) then ! DST-1
                call d_forward_trig_transform(rhs_tt(i,j,:),handle_z,ipar_z,dpar_z,status)
            end if
        end do
    end do

    !if (pbc_z==2) then
    !    rhs_tt(1:nx,1:ny,1:nz)=rhs_tt(1:nx,1:ny,1:nz)/eig_tt(:,:,nz:1:-1)
    !else
    !    rhs_tt(1:nx,1:ny,1:nz)=rhs_tt(1:nx,1:ny,1:nz)/eig_tt
    !end if
    rhs_tt(1:nx,1:ny,1:nz)=rhs_tt(1:nx,1:ny,1:nz)/eig_tt

    !z-direction
    do j=1,ny
        do i=1,nx
            if (pbc_z==2) then ! iDST-2 = reverse order of input, iDCT-3, flip sign of every other output
                !rhs_tt(i,j,1:nz)=rhs_tt(i,j,nz:1:-1)
                call d_forward_trig_transform(rhs_tt(i,j,:),handle_z,ipar_z,dpar_z,status)
                rhs_tt(i,j,1:nz:2)=-rhs_tt(i,j,1:nz:2)
            elseif (pbc_z==3) then ! iDCT-2 = DCT-3
                call d_forward_trig_transform(rhs_tt(i,j,:),handle_z,ipar_z,dpar_z,status)
            elseif (pbc_z==4) then ! iDST-1
                call d_backward_trig_transform(rhs_tt(i,j,:),handle_z,ipar_z,dpar_z,status)
            end if
        end do
    end do

    !y-direction
    do k=1,nz
        do i=1,nx
            if (pbc_y==2) then ! iDST-2 = reverse order of input, iDCT-3, flip sign of every other output
                rhs_tt(i,1:ny,k)=rhs_tt(i,ny:1:-1,k)
                call d_forward_trig_transform(rhs_tt(i,:,k),handle_y,ipar_y,dpar_y,status)
                rhs_tt(i,1:ny:2,k)=-rhs_tt(i,1:ny:2,k)
            elseif (pbc_y==3) then ! iDCT-2 = DCT-3
                call d_forward_trig_transform(rhs_tt(i,:,k),handle_y,ipar_y,dpar_y,status)
            elseif (pbc_y==4) then ! iDST-1
                call d_backward_trig_transform(rhs_tt(i,:,k),handle_y,ipar_y,dpar_y,status)
            end if
        end do
    end do

    !x-direction
    do k=1,nz
        do j=1,ny
            if (pbc_x==2) then ! iDST-2 = reverse order of input, iDCT-3, flip sign of every other output
                rhs_tt(1:nx,j,k)=rhs_tt(nx:1:-1,j,k)
                call d_forward_trig_transform(rhs_tt(:,j,k),handle_x,ipar_x,dpar_x,status)
                rhs_tt(1:nx:2,j,k)=-rhs_tt(1:nx:2,j,k)
            elseif (pbc_x==3) then ! iDCT-2 = DCT-3
                call d_forward_trig_transform(rhs_tt(:,j,k),handle_x,ipar_x,dpar_x,status)
            elseif (pbc_x==4) then ! iDST-1
                call d_backward_trig_transform(rhs_tt(:,j,k),handle_x,ipar_x,dpar_x,status)
            end if
        end do
    end do

    end subroutine DCT_poisson_solver

    end module NS_functions

