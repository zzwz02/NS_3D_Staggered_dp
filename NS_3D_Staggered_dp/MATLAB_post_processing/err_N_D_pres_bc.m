clear all
close all
clc

color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2]; t5=[1:dt(5):2];
err_1 = import_err_file('../err_file/err_vel_AB2_3.E-4_D_Ustar_TOffset.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_Ustar_TOffset.txt');

pcolor=color(1,:);
figure(1)
semilogy([0:4000],err_1.rel_err_u,'-','Color',color(1,:))
hold on;
plot([0:4000],err_2.rel_err_u,'--','Color',color(1,:))
plot([0:4000],err_1.rel_err_p,'-','Color',color(2,:))
plot([0:4000],err_2.rel_err_p,'--','Color',color(2,:))

legend('u error with D pres. bc','u error with N pres. bc','p error with D pres. bc','p error with N pres. bc')
xlabel('time step'); ylabel('L_\infty rel. error')