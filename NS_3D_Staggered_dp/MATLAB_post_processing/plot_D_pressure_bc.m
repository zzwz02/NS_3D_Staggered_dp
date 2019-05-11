clear all
clc

t_start=0; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13;0.49,0.18,0.56];
dt=[4e-3,4e-3,2e-3,2e-3,1e-3,1e-3,5e-4,5e-4]; 
t1=[t_start:dt(1)*2:t_end]; t2=[t_start:dt(2)*2:t_end]; t3=[t_start:dt(3)*2:t_end]; t4=[t_start:dt(4)*2:t_end];
t5=[t_start:dt(5)*2:t_end]; t6=[t_start:dt(6)*2:t_end]; t7=[t_start:dt(7)*2:t_end]; t8=[t_start:dt(8)*2:t_end];

close all
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_N_U_TOffset_restart.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_D_U_TOffset_restart.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_2.E-3_p_N_U_TOffset_restart.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_2.E-3_p_D_U_TOffset_restart.txt');
err_5 = import_err_file('../err_file/err_vel_AB2_1.E-3_p_N_U_TOffset_restart.txt');
err_6 = import_err_file('../err_file/err_vel_AB2_1.E-3_p_D_U_TOffset_restart.txt');
err_7 = import_err_file('../err_file/err_vel_AB2_5.E-4_p_N_U_TOffset_restart.txt');
err_8 = import_err_file('../err_file/err_vel_AB2_5.E-4_p_D_U_TOffset_restart.txt');

figure
semilogy(t1,err_1.rel_err_u,'-','Color',color(1,:))
hold on;
semilogy(t2,err_2.rel_err_u,'--','Color',color(1,:))
semilogy(t3,err_3.rel_err_u,'-','Color',color(2,:))
semilogy(t4,err_4.rel_err_u,'--','Color',color(2,:))
semilogy(t5,err_5.rel_err_u,'-','Color',color(3,:))
semilogy(t6,err_6.rel_err_u,'--','Color',color(3,:))
semilogy(t7,err_7.rel_err_u,'-','Color',color(4,:))
semilogy(t8,err_8.rel_err_u,'--','Color',color(4,:))

xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{u,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

figure
semilogy(err_1.rel_err_p,'-','Color',color(1,:))
hold on;
semilogy(err_2.rel_err_p,'--','Color',color(1,:))
semilogy(err_3.rel_err_p,'-','Color',color(2,:))
semilogy(err_4.rel_err_p,'--','Color',color(2,:))
semilogy(err_5.rel_err_p,'-','Color',color(3,:))
semilogy(err_6.rel_err_p,'--','Color',color(3,:))
semilogy(err_7.rel_err_p,'-','Color',color(4,:))
semilogy(err_8.rel_err_p,'--','Color',color(4,:))


h1=plot([1 2],[10 10],'-k');
h2=plot([1 2],[10 10],'--k');
legend([h1 h2],'Neumann pressure b.c', 'Dirichlet pressure b.c');
ylim([3e-5 1.2e-1])

xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{p,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

%%
dt=[4e-3,2e-3,1e-3,5e-4]; 
err_1 = import_err_file('../err_file/err_vel_AB2_4.E-3_p_D_U_TOffset_restart.txt');
err_2 = import_err_file('../err_file/err_vel_AB2_2.E-3_p_D_U_TOffset_restart.txt');
err_3 = import_err_file('../err_file/err_vel_AB2_1.E-3_p_D_U_TOffset_restart.txt');
err_4 = import_err_file('../err_file/err_vel_AB2_5.E-4_p_D_U_TOffset_restart.txt');

err_u_max=[max(err_1.rel_err_u(1:end)) max(err_2.rel_err_u(1:end)) max(err_3.rel_err_u(1:end)) max(err_4.rel_err_u(1:end))];
err_v_max=[max(err_1.rel_err_v(1:end)) max(err_2.rel_err_v(1:end)) max(err_3.rel_err_v(1:end)) max(err_4.rel_err_v(1:end))];
err_w_max=[max(err_1.rel_err_w(1:end)) max(err_2.rel_err_w(1:end)) max(err_3.rel_err_w(1:end)) max(err_4.rel_err_w(1:end))];
err_p_max=[max(err_1.rel_err_p(1:end)) max(err_2.rel_err_p(1:end)) max(err_3.rel_err_p(1:end)) max(err_4.rel_err_p(1:end))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^1*9e0,'--','Color',[0.47,0.67,0.19])
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\sigma$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
xlim([1e-4 1e-2])
%legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);

err_u_max=[max(err_1.rel_err_u(t1>=1.02)) max(err_2.rel_err_u(t2>=1.02)) max(err_3.rel_err_u(t3>=1.02)) max(err_4.rel_err_u(t4>=1.02))];
err_v_max=[max(err_1.rel_err_v(t1>=1.02)) max(err_2.rel_err_v(t2>=1.02)) max(err_3.rel_err_v(t3>=1.02)) max(err_4.rel_err_v(t4>=1.02))];
err_w_max=[max(err_1.rel_err_w(t1>=1.02)) max(err_2.rel_err_w(t2>=1.02)) max(err_3.rel_err_w(t3>=1.02)) max(err_4.rel_err_w(t4>=1.02))];
err_p_max=[max(err_1.rel_err_p(t1>=1.02)) max(err_2.rel_err_p(t2>=1.02)) max(err_3.rel_err_p(t3>=1.02)) max(err_4.rel_err_p(t4>=1.02))];

figure
loglog(dt, err_u_max,'^')
hold on;
loglog(dt, err_v_max,'s')
loglog(dt, err_w_max,'o')
loglog(dt, err_p_max,'+')
plot([2e-4 7e-3],[2e-4 7e-3].^2*3e2,'-.','Color',[0.47,0.67,0.19])

xlabel('$\sigma$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
xlim([1e-4 1e-2])
legend('u','v','w','p')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);