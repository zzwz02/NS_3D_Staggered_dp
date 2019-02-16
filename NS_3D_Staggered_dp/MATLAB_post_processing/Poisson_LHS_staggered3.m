function LHS_poisson = Poisson_LHS_staggered3(nxp,nyp,nzp,dx2,dy2,dz2,pbc_x,pbc_y,pbc_z,dx,dy,dz,decomp)

if (nzp<1)
    nzp=1;
end

Ix=speye(nxp);
Iy=speye(nyp);
Iz=speye(nzp);
e = ones(nxp,1); Ax = spdiags([e -2*e e], -1:1, nxp, nxp)/(dx2);
e = ones(nyp,1); Ay = spdiags([e -2*e e], -1:1, nyp, nyp)/(dy2);
e = ones(nzp,1); Az = spdiags([e -2*e e], -1:1, nzp, nzp)/(dz2);
if (pbc_x==1)
    Ax(1,nxp) = 1/(dx2);
    Ax(nxp,1) = 1/(dx2);
elseif (pbc_x==2)
    Ax(1,1) = -3/(dx2);
    Ax(nxp,nxp) = -3/(dx2);
elseif (pbc_x==3)
    Ax(1,1) = -1/(dx2);
    Ax(nxp,nxp) = -1/(dx2);
end
if (pbc_y==1)
    Ay(1,nyp) = 1/(dy2);
    Ay(nyp,1) = 1/(dy2);
elseif (pbc_y==2)
    Ay(1,1) = -3/(dy2);
    Ay(nyp,nyp) = -3/(dy2);
elseif (pbc_y==3)
    Ay(1,1) = -1/(dy2);
    Ay(nyp,nyp) = -1/(dy2);
end
if (pbc_z==1)
    Az(1,nzp) = 1/(dz2);
    Az(nzp,1) = 1/(dz2);
elseif (pbc_z==2)
    Az(1,1) = -3/(dz2);
    Az(nzp,nzp) = -3/(dz2);
elseif (pbc_z==3)
    Az(1,1) = -1/(dz2);
    Az(nzp,nzp) = -1/(dz2);
end

if (nzp==1)
    LHS_poisson = kron(Iy,Ax)+kron(Ay,Ix);
else
    LHS_poisson = kron(kron(Iz,Iy),Ax)+kron(kron(Iz,Ay),Ix)+kron(kron(Az,Iy),Ix);
end

LHS_poisson=LHS_poisson;

% LHS_poisson = kron(kron(Iz,Iy),K1(nx,dx,1))+kron(kron(Iz,K1(ny,dy,1)),Ix)+kron(kron(K1(nz,dz,1),Iy),Ix);
% LHS_poisson(1,1) = 3/2*LHS_poisson(1,1);
%LHS_poisson(1,1) = 3/2*LHS_poisson(1,1);

if (decomp)
    LHS_poisson = decomposition(LHS_poisson,'lu');
end
end


function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
end