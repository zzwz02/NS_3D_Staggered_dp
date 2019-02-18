clear all
%close all
clc

t_start=1; t_end=2;
dt=[4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4];
dt_string={'4.E-3','2.E-3','1.E-3','5.E-4','3.E-4'};

idx=1;
t=[t_start:dt(idx):t_end];
err_dp_max=zeros(length(t),1);
resim_file=['../err_file/err_vel_AB2_',dt_string{idx},'_N_U_TOffset_restart.h5'];
bigDNS_file=['../',h5readatt(resim_file,'/','big_DNS_file')];
slice=16;

temp=h5info(resim_file);

for i=2:length(t)
    gName=temp.Groups(i).Name;
%     u_star=h5read(resim_file,[gName,'/u_star_sub']); %u=u(:,2:end-1,2:end-1);
%     RHS=h5read(resim_file,[gName,'/RHS_poisson_sub']); RHS=RHS(2:end-1,2:end-1,2:end-1);
    dp=h5read(resim_file,[gName,'/dp_sub']); %dp=dp-mean(dp(:));
%     p=h5read(resim_file,[gName,'/p_sub']); %p=p-mean(p(:));
%     u=h5read(resim_file,[gName,'/u_sub']); %u=u(:,2:end-1,2:end-1);
    
%     u_star_sub=h5read(bigDNS_file,[gName,'/u_star_sub']); %u_sub=u_sub(:,2:end-1,2:end-1);
%     RHS_sub=h5read(bigDNS_file,[gName,'/RHS_poisson_sub']); RHS_sub=RHS_sub(2:end-1,2:end-1,2:end-1);
    p_sub=h5read(bigDNS_file,[gName,'/p_sub']); %p_sub=p_sub-mean(p_sub(:));
    dp_sub=h5read(bigDNS_file,[gName,'/dp_sub']); %p_sub=p_sub-mean(p_sub(:));
%     u_sub=h5read(bigDNS_file,[gName,'/u_sub']); %u=u(:,2:end-1,2:end-1);;
    
%     err_u_star=u_star-u_star_sub;
%     err_RHS=RHS-RHS_sub;
    err_dp=dp-dp_sub;
%     err_p=p-p_sub;
%     err_u=u-u_sub;
    
%     err_u_star_max(i)=max(abs(err_u_star(:)));
%     err_RHS_max(i)=max(abs(err_RHS(:)));
    err_dp_max(i)=max(abs(err_dp(:)))/rms(dp_sub(:));
%     err_p_max(i)=max(abs(err_p(:)))/rms(p_sub(:));
%     err_u_max(i)=max(abs(err_u(:)));
    
end

figure(1);
semilogy(t,err_dp_max,'-+')
hold on;
