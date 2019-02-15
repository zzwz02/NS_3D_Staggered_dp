clear all
close all
clc

color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13;0.49,0.18,0.56];

dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; sub_tstep=[1,2,3,4,5,10,20,50,100];
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2]; t5=[1:dt(5):2];
err_1 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset.txt');
err_2 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep2.txt');
err_3 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep3.txt');
err_4 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep4.txt');
err_5 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep5.txt');
err_6 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep10.txt');
err_7 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep20.txt');
err_8 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep50.txt');
err_9 = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep100.txt');

tt=t2;
pcolor=color(1,:);
figure(4)
if (exist('err_1','var')) semilogy(tt,err_1.rel_err_p,'--'); end
hold on;
if (exist('err_2','var')) semilogy(tt,err_2.rel_err_p); end
if (exist('err_5','var')) semilogy(tt,err_5.rel_err_p); end
if (exist('err_6','var')) semilogy(tt,err_6.rel_err_p); end
if (exist('err_8','var')) semilogy(tt,err_8.rel_err_p); end

crit=1;
if (exist('err_1','var')) err_u_max(1)=max(err_1.rel_err_u(11:end)); else err_u_max(1)=NaN; end
if (exist('err_2','var')) err_u_max(2)=max(err_2.rel_err_u(11:end)); else err_u_max(2)=NaN; end
if (exist('err_3','var')) err_u_max(3)=max(err_3.rel_err_u(11:end)); else err_u_max(3)=NaN; end
if (exist('err_4','var')) err_u_max(4)=max(err_4.rel_err_u(11:end)); else err_u_max(4)=NaN; end
if (exist('err_5','var')) err_u_max(5)=max(err_5.rel_err_u(11:end)); else err_u_max(5)=NaN; end
if (exist('err_6','var')) err_u_max(6)=max(err_6.rel_err_u(11:end)); else err_u_max(6)=NaN; end
if (exist('err_7','var')) err_u_max(7)=max(err_7.rel_err_u(11:end)); else err_u_max(7)=NaN; end
if (exist('err_8','var')) err_u_max(8)=max(err_8.rel_err_u(11:end)); else err_u_max(8)=NaN; end
if (exist('err_9','var')) err_u_max(9)=max(err_9.rel_err_u(11:end)); else err_u_max(9)=NaN; end

figure(5)
loglog(sub_tstep,err_u_max,'^','Color',color(1,:))
hold on;
temp1=log(dt); temp2=log(err_u_max);

if (exist('err_1','var')) err_v_max(1)=max(err_1.rel_err_v(11:end)); else err_v_max(1)=NaN; end
if (exist('err_2','var')) err_v_max(2)=max(err_2.rel_err_v(11:end)); else err_v_max(2)=NaN; end
if (exist('err_3','var')) err_v_max(3)=max(err_3.rel_err_v(11:end)); else err_v_max(3)=NaN; end
if (exist('err_4','var')) err_v_max(4)=max(err_4.rel_err_v(11:end)); else err_v_max(4)=NaN; end
if (exist('err_5','var')) err_v_max(5)=max(err_5.rel_err_v(11:end)); else err_v_max(5)=NaN; end
if (exist('err_6','var')) err_v_max(6)=max(err_6.rel_err_v(11:end)); else err_v_max(6)=NaN; end
if (exist('err_7','var')) err_v_max(7)=max(err_7.rel_err_v(11:end)); else err_v_max(7)=NaN; end
if (exist('err_8','var')) err_v_max(8)=max(err_8.rel_err_v(11:end)); else err_v_max(8)=NaN; end
if (exist('err_9','var')) err_v_max(9)=max(err_9.rel_err_v(11:end)); else err_v_max(9)=NaN; end

figure(5)
loglog(sub_tstep,err_v_max,'s','Color',color(2,:))
hold on;
temp1=log(dt); temp2=log(err_v_max);

if (exist('err_1','var')) err_w_max(1)=max(err_1.rel_err_w(11:end)); else err_w_max(1)=NaN; end
if (exist('err_2','var')) err_w_max(2)=max(err_2.rel_err_w(11:end)); else err_w_max(2)=NaN; end
if (exist('err_3','var')) err_w_max(3)=max(err_3.rel_err_w(11:end)); else err_w_max(3)=NaN; end
if (exist('err_4','var')) err_w_max(4)=max(err_4.rel_err_w(11:end)); else err_w_max(4)=NaN; end
if (exist('err_5','var')) err_w_max(5)=max(err_5.rel_err_w(11:end)); else err_w_max(5)=NaN; end
if (exist('err_6','var')) err_w_max(6)=max(err_6.rel_err_w(11:end)); else err_w_max(6)=NaN; end
if (exist('err_7','var')) err_w_max(7)=max(err_7.rel_err_w(11:end)); else err_w_max(7)=NaN; end
if (exist('err_8','var')) err_w_max(8)=max(err_8.rel_err_w(11:end)); else err_w_max(8)=NaN; end
if (exist('err_9','var')) err_w_max(9)=max(err_9.rel_err_w(11:end)); else err_w_max(9)=NaN; end

figure(5)
loglog(sub_tstep,err_w_max,'o','Color',color(3,:))
hold on;
temp1=log(dt); temp2=log(err_w_max);

if (exist('err_1','var')) err_p_max(1)=max(err_1.rel_err_p(11:end)); else err_p_max(1)=NaN; end
if (exist('err_2','var')) err_p_max(2)=max(err_2.rel_err_p(11:end)); else err_p_max(2)=NaN; end
if (exist('err_3','var')) err_p_max(3)=max(err_3.rel_err_p(11:end)); else err_p_max(3)=NaN; end
if (exist('err_4','var')) err_p_max(4)=max(err_4.rel_err_p(11:end)); else err_p_max(4)=NaN; end
if (exist('err_5','var')) err_p_max(5)=max(err_5.rel_err_p(11:end)); else err_p_max(5)=NaN; end
if (exist('err_6','var')) err_p_max(6)=max(err_6.rel_err_p(11:end)); else err_p_max(6)=NaN; end
if (exist('err_7','var')) err_p_max(7)=max(err_7.rel_err_p(11:end)); else err_p_max(7)=NaN; end
if (exist('err_8','var')) err_p_max(8)=max(err_8.rel_err_p(11:end)); else err_p_max(8)=NaN; end
if (exist('err_9','var')) err_p_max(9)=max(err_9.rel_err_p(11:end)); else err_p_max(9)=NaN; end

figure(5)
loglog(sub_tstep,err_p_max,'+','Color',color(4,:))
hold on;
loglog([1*0.9 10*1.1],[1*0.9 10*1.1].^-2*1e-4,'k--')
loglog([1*0.9 10*1.1],[1*0.9 10*1.1].^-2*4.5e-4,'k--')
xlim([0.8, 110])
polyfit(log(sub_tstep(1:4)),log(err_u_max(1:4)),1)

legend('u','v','w','p','slope=2')
xlabel('number of sub-time step, dt_0/dt'); ylabel('L_\infty rel. error (compared with AB2 solut.)')
