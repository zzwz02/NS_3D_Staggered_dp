clear all
close all

root_folder='D:\Documents\source\repos\NS_3D_Staggered_dp\NS_3D_Staggered_dp\';
data_folder=[root_folder,'AB2_result\p-form\'];
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
filename=[data_folder,'HIT_256^3_decay_4.E-3_AB2_p_t_.0000.h5'];
disp(filename);

tic
u=h5read(filename,'/u');
v=h5read(filename,'/v');
w=h5read(filename,'/w');
toc;

u_p=(u(1:end-1,2:end-1,2:end-1)+u(2:end,2:end-1,2:end-1))/2;
v_p=(v(2:end-1,1:end-1,2:end-1)+v(2:end-1,2:end,2:end-1))/2;
w_p=(w(2:end-1,2:end-1,1:end-1)+w(2:end-1,2:end-1,2:end))/2;

[spectrum,k,bin_counter,time,r,dx] = PowerSpec(u_p,v_p,w_p,lx,nx);
[Dissipation,kin_Sp,kin_Ph,kin_E,up] = SpecProp(spectrum,k,nu,u_p,v_p,w_p,nx,dx);
[eta,u_eta,tau]=KolmoScale(nu,Dissipation);

2*pi/nx/eta
urms=rms([u_p(:);v_p(:);w_p(:)]);
lambda=sqrt(15*nu*urms^2/Dissipation);
Re_lambda=up*lambda/nu;

figure(1);
loglog(k, spectrum)
hold on;
%plot(k, 4e-1*k.^(-5/3))

figure(2);
loglog(k*eta,2*nu*k.^2.*spectrum'/(Dissipation*nu)^0.75)
hold on;
%%
filename=[data_folder,'HIT_256^3_decay_4.E-3_AB2_p_t_2.0000.h5'];
disp(filename);

tic
u=h5read(filename,'/u');
v=h5read(filename,'/v');
w=h5read(filename,'/w');
toc;

u_p=(u(1:end-1,2:end-1,2:end-1)+u(2:end,2:end-1,2:end-1))/2;
v_p=(v(2:end-1,1:end-1,2:end-1)+v(2:end-1,2:end,2:end-1))/2;
w_p=(w(2:end-1,2:end-1,1:end-1)+w(2:end-1,2:end-1,2:end))/2;

[spectrum,k,bin_counter,time,r,dx] = PowerSpec(u_p,v_p,w_p,lx,nx);
[Dissipation,kin_Sp,kin_Ph,kin_E,up] = SpecProp(spectrum,k,nu,u_p,v_p,w_p,nx,dx);
[eta1,u_eta,tau]=KolmoScale(nu,Dissipation);

2*pi/nx/eta1

figure(1);
loglog(k, spectrum)
hold on;
plot(k, 4e-1*k.^(-5/3))
legend('time step=0.       dx/\eta=1.33','time step=500. dx/\eta=0.58')
xlabel('$k$','Interpreter','latex'); ylabel('$E(k)$','Interpreter','latex')
xlim([0.8 2e2])

figure(2);
loglog(k*eta1,2*nu*k.^2.*spectrum'/(Dissipation*nu)^0.75)
legend('time step=0.       dx/\eta=1.33','time step=500. dx/\eta=0.58')
xlabel('$k\eta$','Interpreter','latex'); ylabel('$k^2 E(k)/(\varepsilon\nu)^{3/4}$','Interpreter','latex')
xlim([1e-2 3e0])
%%
time=0;
filename=[data_folder,'HIT_256^3_decay_4.E-3_AB2_slice_129.h5'];
group_name=num2str(time,'%0.4f');
if (time==0)
    group_name=group_name(2:end);
end

u_slice=h5read(filename,['/t_',group_name,'/u_slice']);
contourf(u_slice)