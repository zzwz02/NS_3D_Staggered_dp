clear all
close all
clc

dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4];
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2]; t5=[1:dt(5):2];
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_D_Ustar_TOffset.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_D_Ustar_TOffset.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_D_Ustar_TOffset.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_D_Ustar_TOffset.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_D_Ustar_TOffset.txt');

figure(1)
semilogy(t1,err_1.rel_err_u)
hold on;
plot(t2,err_2.rel_err_u)
plot(t3,err_3.rel_err_u)
plot(t4,err_4.rel_err_u)
plot(t5,err_5.rel_err_u)

err_u_max(1)=max(err_1.rel_err_u);
err_u_max(2)=max(err_2.rel_err_u);
err_u_max(3)=max(err_3.rel_err_u);
err_u_max(4)=max(err_4.rel_err_u);
err_u_max(5)=max(err_5.rel_err_u);

figure(2)
loglog(dt,err_u_max,'o')
hold on;
temp1=log(dt); temp2=log(err_u_max);