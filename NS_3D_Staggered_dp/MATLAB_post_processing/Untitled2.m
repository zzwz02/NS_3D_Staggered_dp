clear all
%close all
%clc

t_start=1; t_end=2;
dt=[4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4];
dt_string={'4.E-3','2.E-3','1.E-3','5.E-4','3.E-4'};

err_u_star_max=zeros(length(length(dt)),1);
err_RHS_max=err_u_star_max; err_dp_max=err_u_star_max; err_p_max=err_u_star_max; err_u_max=err_u_star_max;

t=[t_start:dt(1):t_end];

bigDNS_file=['../AB2_result.MARCC/HIT_256^3_decay_4.E-3_AB2_dp_x0_16_nx0_32_sub.h5'];
%temp=h5info(bigDNS_file);
gName='/t_.0040';


u_star_sub=h5read(bigDNS_file,[gName,'/u_star_sub']); %u_sub=u_sub(:,2:end-1,2:end-1);
RHS_sub=h5read(bigDNS_file,[gName,'/RHS_poisson_sub']); RHS_sub=RHS_sub(2:end-1,2:end-1,2:end-1);
dp_sub=h5read(bigDNS_file,[gName,'/dp_sub']); %dp_sub=dp_sub-mean(dp_sub(:));
p_sub=h5read(bigDNS_file,[gName,'/p_sub']); %p_sub=p_sub-mean(p_sub(:));
u_sub=h5read(bigDNS_file,[gName,'/u_sub']); %u=u(:,2:end-1,2:end-1);

