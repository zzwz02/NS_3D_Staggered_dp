clear all
close all
clc

color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4]; 
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2];% t5=[1:dt(5):2];
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_Ustar_TOffset.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_p_N_Ustar_TOffset.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_p_N_Ustar_TOffset.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_p_N_Ustar_TOffset.txt');
%err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_p_N_Ustar_TOffset.txt');

pcolor=color(1,:);
figure(1)
semilogy(t1,err_1.rel_err_dp,'--')
hold on;
plot(t2,err_2.rel_err_dp,'--')
plot(t3,err_3.rel_err_dp,'--')
plot(t4,err_4.rel_err_dp,'--')
% plot(t5,err_5.rel_err_dp,'--')

%dt=[1e-1 1e-3 1e-5 1e-7 1e-9];
crit=1.004;
err_u_max(1)=max(err_1.rel_err_u(6:end));
err_u_max(2)=max(err_2.rel_err_u(6:end));
err_u_max(3)=max(err_3.rel_err_u(6:end));
err_u_max(4)=max(err_4.rel_err_u(6:end));
% err_u_max(5)=max(err_5.rel_err_u(6:end));

figure(2)
loglog(dt,err_u_max,'^','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_v(6:end));
err_u_max(2)=max(err_2.rel_err_v(6:end));
err_u_max(3)=max(err_3.rel_err_v(6:end));
err_u_max(4)=max(err_4.rel_err_v(6:end));
% err_u_max(5)=max(err_5.rel_err_v(6:end));

figure(2)
loglog(dt,err_u_max,'s','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_w(6:end));
err_u_max(2)=max(err_2.rel_err_w(6:end));
err_u_max(3)=max(err_3.rel_err_w(6:end));
err_u_max(4)=max(err_4.rel_err_w(6:end));
% err_u_max(5)=max(err_5.rel_err_w(6:end));

figure(2)
loglog(dt,err_u_max,'o','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_dp(6:end));
err_u_max(2)=max(err_2.rel_err_dp(6:end));
err_u_max(3)=max(err_3.rel_err_dp(6:end));
err_u_max(4)=max(err_4.rel_err_dp(6:end));
% err_u_max(5)=max(err_5.rel_err_dp(6:end));

figure(2)
loglog(dt,err_u_max,'+','Color',pcolor)
hold on;
plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^2*3e2,'--','Color',[0.47,0.67,0.19])
plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^1*9e0,'-','Color',[0.47,0.67,0.19])
%%
err_1 = import_err_file('../err_file/err_vel_Euler_4.E-3_p_N_Ustar_TOffset.txt');
err_2 = import_err_file('../err_file/err_vel_Euler_2.E-3_p_N_Ustar_TOffset.txt');
err_3 = import_err_file('../err_file/err_vel_Euler_1.E-3_p_N_Ustar_TOffset.txt');
err_4 = import_err_file('../err_file/err_vel_Euler_5.E-4_p_N_Ustar_TOffset.txt');
% err_5 = import_err_file('../err_file/err_vel_Euler_3.E-4_p_N_Ustar_TOffset.txt');

pcolor=color(2,:);
figure(1)
semilogy(t1,err_1.rel_err_dp,'--')
hold on;
plot(t2,err_2.rel_err_dp,'--')
plot(t3,err_3.rel_err_dp,'--')
plot(t4,err_4.rel_err_dp,'--')
% plot(t5,err_5.rel_err_dp,'--')

%dt=[1e-1 1e-3 1e-5 1e-7 1e-9];
crit=1.004;
err_u_max(1)=max(err_1.rel_err_u(6:end));
err_u_max(2)=max(err_2.rel_err_u(6:end));
err_u_max(3)=max(err_3.rel_err_u(6:end));
err_u_max(4)=max(err_4.rel_err_u(6:end));
% err_u_max(5)=max(err_5.rel_err_u(6:end));

figure(2)
loglog(dt,err_u_max,'^','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_v(6:end));
err_u_max(2)=max(err_2.rel_err_v(6:end));
err_u_max(3)=max(err_3.rel_err_v(6:end));
err_u_max(4)=max(err_4.rel_err_v(6:end));
% err_u_max(5)=max(err_5.rel_err_v(6:end));

figure(2)
loglog(dt,err_u_max,'s','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_w(6:end));
err_u_max(2)=max(err_2.rel_err_w(6:end));
err_u_max(3)=max(err_3.rel_err_w(6:end));
err_u_max(4)=max(err_4.rel_err_w(6:end));
% err_u_max(5)=max(err_5.rel_err_w(6:end));

figure(2)
loglog(dt,err_u_max,'o','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_dp(6:end));
err_u_max(2)=max(err_2.rel_err_dp(6:end));
err_u_max(3)=max(err_3.rel_err_dp(6:end));
err_u_max(4)=max(err_4.rel_err_dp(6:end));
% err_u_max(5)=max(err_5.rel_err_dp(6:end));

figure(2)
loglog(dt,err_u_max,'+','Color',pcolor)
hold on;

plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^2*2.4e2,'--','Color',[0.47,0.67,0.19])
plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^1*3.5e1,'-','Color',[0.47,0.67,0.19])
%plot([dt(1)/0.9 dt(end)*1.1],[dt(1)/0.9 dt(end)*1.1].^-1*1e-6,'--')

annotation('textbox',[0.5,0.5,0.2,0.2],'String',"slope=1, Euler",'FitBoxToText','on','LineStyle','none');
annotation('textbox',[0.2,0.2,0.2,0.2],'String',"slope=2, start with Wuler then AB2",'FitBoxToText','on','LineStyle','none');
xlabel('\Deltat'); ylabel('L_\infty relative error')
polyfit(log(dt),log(err_u_max),1)