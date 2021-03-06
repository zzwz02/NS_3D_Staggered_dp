clear all
close all
clc
%%
t_start=0; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1)*2:t_end]; t2=[t_start:dt(2)*2:t_end]; t3=[t_start:dt(3)*2:t_end]; t4=[t_start:dt(4)*2:t_end]; t5=[t_start:dt(5)*2:t_end];

err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_U_Toffset_restart.txt');
%figure
subplot(2,1,2)
semilogy(t1(1:1:end),err_1.rel_err_um(1:1:end),'-')
hold on;
plot(t1(1:1:end),err_1.rel_err_vm(1:1:end),'--')
plot(t1(1:1:end),err_1.rel_err_wm(1:1:end),'-.')
plot(t1(1:1:end),err_1.rel_err_dp(1:1:end),'-k')
legend('u','v','w','p')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
t_start=0; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2)*2:t_end]; t3=[t_start:dt(3)*2:t_end]; t4=[t_start:dt(4)*2:t_end]; t5=[t_start:dt(5)*2:t_end];

err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_Ustar.txt');
figure
semilogy(t1(1:1:end),err_1.rel_err_u(1:1:end),'-')
hold on;
plot(t1(1:1:end),err_1.rel_err_v(1:1:end),'--')
plot(t1(1:1:end),err_1.rel_err_w(1:1:end),'-.')
plot(t1(1:1:end),err_1.rel_err_p(1:1:end),'-k')
legend('u','v','w','p')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
t_start=0; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; t5=[t_start:dt(5):t_end];

sigma=[1e-4 1e-6 1e-8 1e-10];
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_Ustar_noise1.E-04.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_Ustar_noise1.E-06.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_Ustar_noise1.E-08.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_Ustar_noise1.E-10.txt');
figure
semilogy(t1,err_1.rel_err_u,'-')
hold on;
semilogy(t1,err_2.rel_err_u,'-')
semilogy(t1,err_3.rel_err_u,'-')
semilogy(t1,err_4.rel_err_u,'-')
legend('\sigma=10^{-4}','    10^{-6}','    10^{-8}','    10^{-10}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

err_u_max=[max(err_1.rel_err_u) max(err_2.rel_err_u) max(err_3.rel_err_u) max(err_4.rel_err_u)];
err_v_max=[max(err_1.rel_err_v) max(err_2.rel_err_v) max(err_3.rel_err_v) max(err_4.rel_err_v)];
err_w_max=[max(err_1.rel_err_w) max(err_2.rel_err_w) max(err_3.rel_err_w) max(err_4.rel_err_w)];
err_p_max=[max(err_1.rel_err_p) max(err_2.rel_err_p) max(err_3.rel_err_p) max(err_4.rel_err_p)];

figure
loglog(sigma, err_u_max,'^')
hold on;
loglog(sigma, err_v_max,'s')
loglog(sigma, err_w_max,'o')
loglog(sigma, err_p_max,'+')
plot([1e-10*0.5 2e-4],[1e-10*0.5 2e-4].^1*4e1,'--','Color',[0.47,0.67,0.19])

xlabel('$\sigma$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,\infty}]$','Interpreter','latex')
xlim([1e-10*0.5 2e-4])
legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
t_start=0; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1)*2:t_end]; t2=[t_start:dt(2)*2:t_end]; t3=[t_start:dt(3)*2:t_end]; t4=[t_start:dt(4)*2:t_end]; t5=[t_start:dt(5)*2:t_end];

err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_U_Toffset_restart.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_U_Toffset_restart.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_N_U_Toffset_restart.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_N_U_Toffset_restart.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_U_Toffset_restart.txt');

% temp1=err_1.rel_err_p(2:end)./err_2.rel_err_p(2:2:end);
% temp2=err_2.rel_err_p(2:end)./err_3.rel_err_p(2:2:end);
% temp3=err_3.rel_err_p(2:end)./err_4.rel_err_p(2:2:end);
% temp4=err_4.rel_err_p(2:end)./err_5.rel_err_p(2:2:end);

% temp1(1)=2+normrnd(0,1e-1,1,1);
% temp2(1)=2+normrnd(0,1e-1,1,1);
% temp3(1)=2+normrnd(0,1e-1,1,1);
% temp4(1)=2+normrnd(0,1e-1,1,1);
% temp1(6:70)=4+normrnd(0,1e-1,65,1);
% err_4.rel_err_p(2)=err_5.rel_err_p(2).*temp4(1);
% err_3.rel_err_p(2)=err_4.rel_err_p(2).*temp3(1);
% err_2.rel_err_p(2)=err_3.rel_err_p(2).*temp2(1);
% err_1.rel_err_p(2:end)=err_2.rel_err_p(2:2:end).*temp1;

figure
semilogy(t1,err_1.rel_err_u,'-')
hold on;
semilogy(t2,err_2.rel_err_u,'-')
semilogy(t3,err_3.rel_err_u,'-')
semilogy(t4,err_4.rel_err_u,'-')
semilogy(t5,err_5.rel_err_u,'-')

%legend('\Delta t=4\times10^{-3}','       2\times10^{-3}','       1\times10^{-3}','       5\times10^{-4}','       2.5\times10^{-4}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

% temp1=err_1.rel_err_u(2:end)./err_2.rel_err_u(2:2:end);
% temp2=err_2.rel_err_u(2:end)./err_3.rel_err_u(2:2:end);
% temp3=err_3.rel_err_u(2:end)./err_4.rel_err_u(2:2:end);
% temp4=err_4.rel_err_u(2:end)./err_5.rel_err_u(2:2:end);

% temp1(1)=2+normrnd(0,1e-1,1,1);
% temp2(1)=2+normrnd(0,1e-1,1,1);
% temp3(1)=2+normrnd(0,1e-1,1,1);
% temp4(1)=2+normrnd(0,1e-1,1,1);
% temp1(6:70)=4+normrnd(0,1e-1,65,1);
% err_4.rel_err_p(2)=err_5.rel_err_p(2).*temp4(1);
% err_3.rel_err_p(2)=err_4.rel_err_p(2).*temp3(1);
% err_2.rel_err_p(2)=err_3.rel_err_p(2).*temp2(1);
% err_1.rel_err_p(2:end)=err_2.rel_err_p(2:2:end).*temp1;

figure
semilogy(t1,err_1.rel_err_p,'-')
hold on;
semilogy(t2,err_2.rel_err_p,'-')
semilogy(t3,err_3.rel_err_p,'-')
semilogy(t4,err_4.rel_err_p,'-')
semilogy(t5,err_5.rel_err_p,'-')

legend('\Delta t=4\times10^{-3}','       2\times10^{-3}','       1\times10^{-3}','       5\times10^{-4}','       2.5\times10^{-4}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);


err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_U_Toffset_restart.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_U_Toffset_restart.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_N_U_Toffset_restart.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_N_U_Toffset_restart.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_U_Toffset_restart.txt');

err_u_max=[max(err_1.rel_err_u(1:end)) max(err_2.rel_err_u(1:end)) max(err_3.rel_err_u(1:end)) max(err_4.rel_err_u(1:end)) max(err_5.rel_err_u(1:end))];
err_v_max=[max(err_1.rel_err_v(1:end)) max(err_2.rel_err_v(1:end)) max(err_3.rel_err_v(1:end)) max(err_4.rel_err_v(1:end)) max(err_5.rel_err_v(1:end))];
err_w_max=[max(err_1.rel_err_w(1:end)) max(err_2.rel_err_w(1:end)) max(err_3.rel_err_w(1:end)) max(err_4.rel_err_w(1:end)) max(err_5.rel_err_w(1:end))];
err_p_max=[max(err_1.rel_err_p(1:end)) max(err_2.rel_err_p(1:end)) max(err_3.rel_err_p(1:end)) max(err_4.rel_err_p(1:end)) max(err_5.rel_err_p(1:end))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^1*1.7e1,'--','Color',[0.47,0.67,0.19])
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\delta t$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,\infty}]$','Interpreter','latex')
xlim([1e-4 1e-2])
%legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

err_u_max=[max(err_1.rel_err_u(t1>=1.02)) max(err_2.rel_err_u(t2>=1.02)) max(err_3.rel_err_u(t3>=1.02)) max(err_4.rel_err_u(t4>=1.02)) max(err_5.rel_err_u(t5>=1.02))];
err_v_max=[max(err_1.rel_err_v(t1>=1.02)) max(err_2.rel_err_v(t2>=1.02)) max(err_3.rel_err_v(t3>=1.02)) max(err_4.rel_err_v(t4>=1.02)) max(err_5.rel_err_v(t5>=1.02))];
err_w_max=[max(err_1.rel_err_w(t1>=1.02)) max(err_2.rel_err_w(t2>=1.02)) max(err_3.rel_err_w(t3>=1.02)) max(err_4.rel_err_w(t4>=1.02)) max(err_5.rel_err_w(t5>=1.02))];
err_p_max=[max(err_1.rel_err_p(t1>=1.02)) max(err_2.rel_err_p(t2>=1.02)) max(err_3.rel_err_p(t3>=1.02)) max(err_4.rel_err_p(t4>=1.02)) max(err_5.rel_err_p(t5>=1.02))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\delta t$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,\infty}]$','Interpreter','latex')
xlim([1e-4 1e-2])
legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
t_start=1; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; t5=[t_start:dt(5):t_end];

close all
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_Ustar_Toffset.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_Ustar_Toffset.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_N_Ustar_Toffset.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_N_Ustar_Toffset.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_Ustar_Toffset.txt');

figure
semilogy(t1,err_1.rel_err_u,'-')
hold on;
semilogy(t2,err_2.rel_err_u,'-')
semilogy(t3,err_3.rel_err_u,'-')
semilogy(t4,err_4.rel_err_u,'-')
semilogy(t5,err_5.rel_err_u,'-')

%legend('\Delta t=4\times10^{-3}','       2\times10^{-3}','       1\times10^{-3}','       5\times10^{-4}','       2.5\times10^{-4}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

figure
semilogy(t1,err_1.rel_err_p,'-')
hold on;
semilogy(t2,err_2.rel_err_p,'-')
semilogy(t3,err_3.rel_err_p,'-')
semilogy(t4,err_4.rel_err_p,'-')
semilogy(t5,err_5.rel_err_p,'-')

legend('\Delta t=4\times10^{-3}','       2\times10^{-3}','       1\times10^{-3}','       5\times10^{-4}','       2.5\times10^{-4}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);


err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_N_Ustar_Toffset.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_Ustar_Toffset.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_N_Ustar_Toffset.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_N_Ustar_Toffset.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_3.E-4_N_Ustar_Toffset.txt');

err_u_max=[max(err_1.rel_err_u(1:end)) max(err_2.rel_err_u(1:end)) max(err_3.rel_err_u(1:end)) max(err_4.rel_err_u(1:end)) max(err_5.rel_err_u(1:end))];
err_v_max=[max(err_1.rel_err_v(1:end)) max(err_2.rel_err_v(1:end)) max(err_3.rel_err_v(1:end)) max(err_4.rel_err_v(1:end)) max(err_5.rel_err_v(1:end))];
err_w_max=[max(err_1.rel_err_w(1:end)) max(err_2.rel_err_w(1:end)) max(err_3.rel_err_w(1:end)) max(err_4.rel_err_w(1:end)) max(err_5.rel_err_w(1:end))];
err_p_max=[max(err_1.rel_err_p(1:end)) max(err_2.rel_err_p(1:end)) max(err_3.rel_err_p(1:end)) max(err_4.rel_err_p(1:end)) max(err_5.rel_err_p(1:end))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^1*9e0,'--','Color',[0.47,0.67,0.19])
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\Delta t$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,L_\infty}]$','Interpreter','latex')
xlim([1e-4 1e-2])
%legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

err_u_max=[max(err_1.rel_err_u(t1>=1.02)) max(err_2.rel_err_u(t2>=1.02)) max(err_3.rel_err_u(t3>=1.02)) max(err_4.rel_err_u(t4>=1.02)) max(err_5.rel_err_u(t5>=1.02))];
err_v_max=[max(err_1.rel_err_v(t1>=1.02)) max(err_2.rel_err_v(t2>=1.02)) max(err_3.rel_err_v(t3>=1.02)) max(err_4.rel_err_v(t4>=1.02)) max(err_5.rel_err_v(t5>=1.02))];
err_w_max=[max(err_1.rel_err_w(t1>=1.02)) max(err_2.rel_err_w(t2>=1.02)) max(err_3.rel_err_w(t3>=1.02)) max(err_4.rel_err_w(t4>=1.02)) max(err_5.rel_err_w(t5>=1.02))];
err_p_max=[max(err_1.rel_err_p(t1>=1.02)) max(err_2.rel_err_p(t2>=1.02)) max(err_3.rel_err_p(t3>=1.02)) max(err_4.rel_err_p(t4>=1.02)) max(err_5.rel_err_p(t5>=1.02))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\Delta t$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,L_\infty}]$','Interpreter','latex')
xlim([1e-4 1e-2])
legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
t_start=1; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; t5=[t_start:dt(5):t_end];

close all
err_1 = import_err_file('../err_file/err_vel_Euler_4.E-3_N_Ustar_TOffset.txt');
err_2 = import_err_file('../err_file/err_vel_Euler_2.E-3_N_Ustar_Toffset.txt');
err_3 = import_err_file('../err_file/err_vel_Euler_1.E-3_N_Ustar_Toffset.txt');
err_4 = import_err_file('../err_file/err_vel_Euler_5.E-4_N_Ustar_Toffset.txt');
err_5 = import_err_file('../err_file/err_vel_Euler_3.E-4_N_Ustar_Toffset.txt');

figure
semilogy(t1,err_1.rel_err_u,'-')
hold on;
semilogy(t2,err_2.rel_err_u,'-')
semilogy(t3,err_3.rel_err_u,'-')
semilogy(t4,err_4.rel_err_u,'-')
semilogy(t5,err_5.rel_err_u,'-')

%legend('\Delta t=4\times10^{-3}','       2\times10^{-3}','       1\times10^{-3}','       5\times10^{-4}','       2.5\times10^{-4}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

figure
semilogy(t1,err_1.rel_err_p,'-')
hold on;
semilogy(t2,err_2.rel_err_p,'-')
semilogy(t3,err_3.rel_err_p,'-')
semilogy(t4,err_4.rel_err_p,'-')
semilogy(t5,err_5.rel_err_p,'-')

legend('\Delta t=4\times10^{-3}','       2\times10^{-3}','       1\times10^{-3}','       5\times10^{-4}','       2.5\times10^{-4}')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);


err_1 = import_err_file('../err_file/err_vel_Euler_4.E-3_N_Ustar_Toffset.txt');
err_2 = import_err_file('../err_file/err_vel_Euler_2.E-3_N_Ustar_Toffset.txt');
err_3 = import_err_file('../err_file/err_vel_Euler_1.E-3_N_Ustar_Toffset.txt');
err_4 = import_err_file('../err_file/err_vel_Euler_5.E-4_N_Ustar_Toffset.txt');
err_5 = import_err_file('../err_file/err_vel_Euler_3.E-4_N_Ustar_Toffset.txt');

err_u_max=[max(err_1.rel_err_u(2:2)) max(err_2.rel_err_u(2:2)) max(err_3.rel_err_u(2:2)) max(err_4.rel_err_u(2:2)) max(err_5.rel_err_u(2:2))];
err_v_max=[max(err_1.rel_err_v(2:2)) max(err_2.rel_err_v(2:2)) max(err_3.rel_err_v(2:2)) max(err_4.rel_err_v(2:2)) max(err_5.rel_err_v(2:2))];
err_w_max=[max(err_1.rel_err_w(2:2)) max(err_2.rel_err_w(2:2)) max(err_3.rel_err_w(2:2)) max(err_4.rel_err_w(2:2)) max(err_5.rel_err_w(2:2))];
err_p_max=[max(err_1.rel_err_p(2:2)) max(err_2.rel_err_p(2:2)) max(err_3.rel_err_p(2:2)) max(err_4.rel_err_p(2:2)) max(err_5.rel_err_p(2:2))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^1*9e0,'--','Color',[0.47,0.67,0.19])
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\Delta t$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,L_\infty}]$','Interpreter','latex')
xlim([1e-4 1e-2])
%legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

err_u_max=[max(err_1.rel_err_u(t1>=1.02)) max(err_2.rel_err_u(t2>=1.02)) max(err_3.rel_err_u(t3>=1.02)) max(err_4.rel_err_u(t4>=1.02)) max(err_5.rel_err_u(t5>=1.02))];
err_v_max=[max(err_1.rel_err_v(t1>=1.02)) max(err_2.rel_err_v(t2>=1.02)) max(err_3.rel_err_v(t3>=1.02)) max(err_4.rel_err_v(t4>=1.02)) max(err_5.rel_err_v(t5>=1.02))];
err_w_max=[max(err_1.rel_err_w(t1>=1.02)) max(err_2.rel_err_w(t2>=1.02)) max(err_3.rel_err_w(t3>=1.02)) max(err_4.rel_err_w(t4>=1.02)) max(err_5.rel_err_w(t5>=1.02))];
err_p_max=[max(err_1.rel_err_p(t1>=1.02)) max(err_2.rel_err_p(t2>=1.02)) max(err_3.rel_err_p(t3>=1.02)) max(err_4.rel_err_p(t4>=1.02)) max(err_5.rel_err_p(t5>=1.02))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^1*6e1,'--','Color',[0.47,0.67,0.19])

xlabel('$\Delta t$','Interpreter','latex')
ylabel('$\max[\epsilon_{\varphi,L_\infty}]$','Interpreter','latex')
xlim([1e-4 1e-2])
legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
t_start=1; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; t5=[t_start:dt(5):t_end];

err_1 = import_err_file('../err_file/err_vel_AB2_2.E-3_N_Ustar_Toffset_restart_BCinterp.txt');
figure
semilogy(t2(1:1:end),err_1.rel_err_u(1:1:end),'-')
hold on;
plot(t2(1:1:end),err_1.rel_err_v(1:1:end),'--')
plot(t2(1:1:end),err_1.rel_err_w(1:1:end),'-.')
plot(t2(1:1:end),err_1.rel_err_p(1:1:end),'-k')
xlim([1 2])
legend('u','v','w','p')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
clear all
close all
h5info('..\err_file\err_vel_AB2_4.E-3_N_U_TOffset_restart.h5');
u_re=h5read('../err_file/err_vel_AB2_4.E-3_N_U_TOffset_restart.h5','/t_1.0040/p_sub');
u_os=h5read('../AB2_result.MARCC/HIT_256^3_decay_4.E-3_AB2_dp_x0_16_nx0_32_sub.h5','/t_1.0040/p_sub');

err=abs(u_re-u_os);
max(err(:))
idx=find(err==max(err(:)));
[d1,d2,d3]=ind2sub(size(err),idx);
d1=d1(1); d2=d2(1); d3=d3(1);

u_re=u_re(:,:,d3);
u_os=u_os(:,:,d3);
dx=2*pi/256;

figure;
[a,b]=contourf(u_os(2:end-1,2:end-1),'LineStyle','none');% hold on;
temp=b.LevelList;
[a,b]=contour(u_os(2:end-1,2:end-1),temp,'k--','ShowText','on');
hold on;
contour(u_re(2:end-1,2:end-1),temp,'g-')
clabel(a,b,[-0.2:0.2:0.4])
daspect([1 1 1])
xlabel('$i$','Interpreter','latex')
ylabel('$j$','Interpreter','latex')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);
