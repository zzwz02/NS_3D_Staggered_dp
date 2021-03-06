clear all
%close all
clc

t_start=1; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; %t5=[t_start:dt(5):t_end];
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_U_Toffset_Restart.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_U_Toffset_Restart.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_N_U_Toffset_Restart.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_N_U_Toffset_Restart.txt');
%err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_U_Toffset_Restart.txt');

pcolor=color(1,:);
figure(1)
semilogy(t1,err_1.rel_err_p,'--')
hold on;
plot(t2,err_2.rel_err_p,'--')
plot(t3,err_3.rel_err_p,'--')
plot(t4,err_4.rel_err_p,'--')
%plot(t5,err_5.rel_err_p,'--')

%dt=[1e-1 1e-3 1e-5 1e-7 1e-9];
crit=1.004;
err_u_mean(1)=mean(err_1.rel_err_u(10:end));
err_u_mean(2)=mean(err_2.rel_err_u(10:end));
err_u_mean(3)=mean(err_3.rel_err_u(10:end));
err_u_mean(4)=mean(err_4.rel_err_u(10:end));
%err_u_mean(5)=mean(err_5.rel_err_u(10:end));
polyfit(log(dt),log(err_u_mean),1)

figure(2)
loglog(dt,err_u_mean,'^','Color',pcolor)
hold on;

err_u_mean(1)=mean(err_1.rel_err_v(10:end));
err_u_mean(2)=mean(err_2.rel_err_v(10:end));
err_u_mean(3)=mean(err_3.rel_err_v(10:end));
err_u_mean(4)=mean(err_4.rel_err_v(10:end));
%err_u_mean(5)=mean(err_5.rel_err_v(10:end));
polyfit(log(dt),log(err_u_mean),1)

figure(2)
loglog(dt,err_u_mean,'s','Color',pcolor)
hold on;

err_u_mean(1)=mean(err_1.rel_err_w(10:end));
err_u_mean(2)=mean(err_2.rel_err_w(10:end));
err_u_mean(3)=mean(err_3.rel_err_w(10:end));
err_u_mean(4)=mean(err_4.rel_err_w(10:end));
% err_u_mean(5)=mean(err_5.rel_err_w(10:end));
polyfit(log(dt),log(err_u_mean),1)

figure(2)
loglog(dt,err_u_mean,'o','Color',pcolor)
hold on;

err_u_mean(1)=mean(err_1.rel_err_p(10:end));
err_u_mean(2)=mean(err_2.rel_err_p(10:end));
err_u_mean(3)=mean(err_3.rel_err_p(10:end));
err_u_mean(4)=mean(err_4.rel_err_p(10:end));
% err_u_mean(5)=mean(err_5.rel_err_p(10:end));
polyfit(log(dt),log(err_u_mean),1)

figure(2)
loglog(dt,err_u_mean,'+','Color',pcolor)
hold on;
temp1=log10(dt); temp2=log10(err_u_mean);

plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^2*2.2e2,'--','Color',[0.47,0.67,0.19])
plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^1*1.9e1,'--','Color',[0.47,0.67,0.19])

legend('u','v','w','p','slope=2','slope=1')
xlabel('\Deltat'); ylabel('L_\infty relative error')
