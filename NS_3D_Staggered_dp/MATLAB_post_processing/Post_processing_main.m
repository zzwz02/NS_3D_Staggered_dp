clear all
close all

nu=0.002;
lx = 2*pi; ly = 2*pi; lz = 2*pi;
nx=256; ny=nx; nz=nx;
dx=lx/nx; dy=ly/ny; dz=lz/nz;
dx2=dx*dx; dy2=dy*dy; dz2=dz*dz;
xp=[-0.5:nx+0.5]*dx; yp=[-0.5:ny+0.5]*dy; zp=[-0.5:nz+0.5]*dz;
xu=[0:nx]*dx;        yu=[-0.5:ny+0.5]*dy; zu=[-0.5:nz+0.5]*dz;
xv=[-0.5:nx+0.5]*dx; yv=[0:ny]*dy;        zv=[-0.5:nz+0.5]*dz;
xw=[-0.5:nx+0.5]*dx; yw=[-0.5:ny+0.5]*dy; zw=[0:nz]*dz;

bc_x=1; bc_y=1; bc_z=1;
pbc_x=1; pbc_y=1; pbc_z=1;

%%
filename='D:\Documents\source\repos\NS_3D_Staggered_dp\NS_3D_Staggered_dp\AB2_result\HIT_256^3_decay_4.E-3_AB2_dp_t_.0000.dat';
disp(filename);

tic
fileID = fopen(filename);
u = fread(fileID,(nx+1)*(ny+2)*(nz+2),'*double'); u=reshape(u,[nx+1, ny+2, nz+2]);
v = fread(fileID,(nx+2)*(ny+1)*(nz+2),'*double'); v=reshape(v,[nx+2, ny+1, nz+2]);
w = fread(fileID,(nx+2)*(ny+2)*(nz+1),'*double'); w=reshape(w,[nx+2, ny+2, nz+1]);
fclose(fileID);
toc;

u_p=(u(1:end-1,2:end-1,2:end-1)+u(2:end,2:end-1,2:end-1))/2;
v_p=(v(2:end-1,1:end-1,2:end-1)+v(2:end-1,2:end,2:end-1))/2;
w_p=(w(2:end-1,2:end-1,1:end-1)+w(2:end-1,2:end-1,2:end))/2;

[spectrum,k,bin_counter,time,r,dx] = PowerSpec(u_p,v_p,w_p,lx,nx);
[Dissipation,kin_Sp,kin_Ph,kin_E,up] = SpecProp(spectrum,k,nu,u_p,v_p,w_p,nx,dx);
[eta,u_eta,tau]=KolmoScale(nu,Dissipation);

2*pi/nx/eta

subplot(1,2,1)
loglog(k, spectrum)
hold on;
%plot(k, 4e-1*k.^(-5/3))

subplot(1,2,2)
loglog(k*eta,2*nu*k.^2.*spectrum')
hold on;
%%
filename='D:\Documents\source\repos\NS_3D_Staggered_dp\NS_3D_Staggered_dp\AB2_result\HIT_256^3_decay_4.E-3_AB2_dp_t_20.0000.dat';
disp(filename);

tic
fileID = fopen(filename);
u = fread(fileID,(nx+1)*(ny+2)*(nz+2),'*double'); u=reshape(u,[nx+1, ny+2, nz+2]);
v = fread(fileID,(nx+2)*(ny+1)*(nz+2),'*double'); v=reshape(v,[nx+2, ny+1, nz+2]);
w = fread(fileID,(nx+2)*(ny+2)*(nz+1),'*double'); w=reshape(w,[nx+2, ny+2, nz+1]);
fclose(fileID);
toc;

u_p=(u(1:end-1,2:end-1,2:end-1)+u(2:end,2:end-1,2:end-1))/2;
v_p=(v(2:end-1,1:end-1,2:end-1)+v(2:end-1,2:end,2:end-1))/2;
w_p=(w(2:end-1,2:end-1,1:end-1)+w(2:end-1,2:end-1,2:end))/2;

[spectrum,k,bin_counter,time,r,dx] = PowerSpec(u_p,v_p,w_p,lx,nx);
[Dissipation,kin_Sp,kin_Ph,kin_E,up] = SpecProp(spectrum,k,nu,u_p,v_p,w_p,nx,dx);
[eta1,u_eta,tau]=KolmoScale(nu,Dissipation);

2*pi/nx/eta1

subplot(1,2,1)
loglog(k, spectrum)
hold on;
plot(k, 4e-1*k.^(-5/3))
legend('time step=0.       dx/\eta=1.33','time step=5000. dx/\eta=0.58')
xlabel('k'); ylabel('E(k)')

subplot(1,2,2)
loglog(k*eta1,2*nu*k.^2.*spectrum')
legend('time step=0.       dx/\eta=1.33','time step=5000. dx/\eta=0.58')
xlabel('k\eta'); ylabel('D(k)')

