clear all
close all
clc

color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2]; t5=[1:dt(5):2];
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_U_TOffset_restart.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_U_TOffset_restart.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_N_U_TOffset_restart.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_N_U_TOffset_restart.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_U_TOffset_restart.txt');

pcolor=color(1,:);
figure(1)
semilogy(t1,err_1.rel_err_p,'--')
hold on;
plot(t2,err_2.rel_err_p,'--')
plot(t3,err_3.rel_err_p,'--')
plot(t4,err_4.rel_err_p,'--')
plot(t5,err_5.rel_err_p,'--')

%dt=[1e-1 1e-3 1e-5 1e-7 1e-9];
crit=1.004;
err_u_max(1)=max(err_1.rel_err_u(3:end));
err_u_max(2)=max(err_2.rel_err_u(3:end));
err_u_max(3)=max(err_3.rel_err_u(3:end));
err_u_max(4)=max(err_4.rel_err_u(3:end));
err_u_max(5)=max(err_5.rel_err_u(3:end));

figure(2)
loglog(dt,err_u_max,'^','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_v(3:end));
err_u_max(2)=max(err_2.rel_err_v(3:end));
err_u_max(3)=max(err_3.rel_err_v(3:end));
err_u_max(4)=max(err_4.rel_err_v(3:end));
err_u_max(5)=max(err_5.rel_err_v(3:end));

figure(2)
loglog(dt,err_u_max,'s','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_w(3:end));
err_u_max(2)=max(err_2.rel_err_w(3:end));
err_u_max(3)=max(err_3.rel_err_w(3:end));
err_u_max(4)=max(err_4.rel_err_w(3:end));
err_u_max(5)=max(err_5.rel_err_w(3:end));

figure(2)
loglog(dt,err_u_max,'o','Color',pcolor)
hold on;

err_u_max(1)=max(err_1.rel_err_p(3:end));
err_u_max(2)=max(err_2.rel_err_p(3:end));
err_u_max(3)=max(err_3.rel_err_p(3:end));
err_u_max(4)=max(err_4.rel_err_p(3:end));
err_u_max(5)=max(err_5.rel_err_p(3:end));

figure(2)
loglog(dt,err_u_max,'+','Color',pcolor)
hold on;

plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^2*2.2e2,'--','Color',[0.47,0.67,0.19])

legend('u','v','w','p','slope=2')
xlabel('\Deltat'); ylabel('L_\infty relative error')
polyfit(log(dt),log(err_u_max),1)