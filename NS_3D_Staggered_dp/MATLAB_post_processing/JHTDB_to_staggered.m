clear all;
filename='D:\Downloads\iso1024_512_cOrder_Velocity.bin';
fileID = fopen(filename,'r');
ini_vel = fread(fileID,'*single'); ini_vel=reshape(ini_vel,[3,512,512,512]);
fclose(fileID);

filename='D:\Downloads\iso1024_512_cOrder_Pressure.bin';
fileID = fopen(filename,'r');
ini_pr = fread(fileID,'*single'); ini_pr=reshape(ini_pr,[512,512,512]);
fclose(fileID);

%%
nx_file=512; nx=256; ny=nx; nz=nx;

temp01=nx_file/nx; temp02=nx_file/ny; temp03=nx_file/nz;
u(1:nx,2:ny+1,2:nz+1)=ini_vel(1,1:temp01:nx_file,2:temp02:nx_file,2:temp03:nx_file);
u(nx+1,:,:)=u(1,:,:);
u(:,1,:)=u(:,ny+1,:); u(:,ny+2,:)=u(:,2,:);
u(:,:,1)=u(:,:,nz+1); u(:,:,nz+2)=u(:,:,2);

v(2:nx+1,1:ny,2:nz+1)=ini_vel(2,2:temp01:nx_file,1:temp02:nx_file,2:temp03:nx_file);
v(1,:,:)=v(nx+1,:,:); v(ny+2,:,:)=v(2,:,:);
v(:,nx+1,:)=v(:,1,:);
v(:,:,1)=v(:,:,nz+1); v(:,:,nz+2)=v(:,:,2);

w(2:nx+1,2:ny+1,1:nz)=ini_vel(3,2:temp01:nx_file,2:temp02:nx_file,1:temp03:nx_file);
w(1,:,:)=w(nx+1,:,:); w(ny+2,:,:)=w(2,:,:);
w(:,1,:)=w(:,ny+1,:); w(:,ny+2,:)=w(:,2,:);
w(:,:,nz+1)=w(:,:,1);

p(2:nx+1,2:ny+1,2:nz+1)=ini_pr(2:temp01:nx_file,2:temp02:nx_file,2:temp03:nx_file);
p(1,:,:)=p(nx+1,:,:); p(ny+2,:,:)=p(2,:,:);
p(:,1,:)=p(:,ny+1,:); p(:,ny+2,:)=p(:,2,:);
p(:,:,1)=p(:,:,nz+1); p(:,:,nz+2)=p(:,:,2);

filename='D:\Documents\source\repos\NS_3D_Staggered_dp\NS_3D_Staggered_dp\init_from_JHTDB.dat';
fileID = fopen(filename,'w');
fwrite(fileID,[u(:);v(:);w(:);p(:)],'*single');
fclose(fileID);