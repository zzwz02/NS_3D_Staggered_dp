%%
t_start=1; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; t5=[t_start:dt(5):t_end];

err_1 = import_err_file('../err_file/spline/err_vel_AB2_5.E-4_N_Ustar_TOffset_sub_tstep10_BCinterp.txt');
figure
semilogy(t4(1:4:end),err_1.rel_err_u(1:4:end),'-')
hold on;
%plot(t4(1:1:end),err_1.rel_err_v(1:1:end),'--')
%plot(t4(1:1:end),err_1.rel_err_w(1:1:end),'-.')
plot(t4(1:4:end),err_1.rel_err_dp(1:4:end),'-k')
xlim([1 2])
legend('u','v','w','p')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')
font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);
%%
clear all
%close all
h5info('..\err_file\err_vel_AB2_4.E-3_N_U_TOffset_restart.h5');
% oz_re=h5read('../err_file/spline/5.E-4_sub_tstep10_BC.h5','/t_2.0000/oz');
% oz_os=h5read('../err_file/spline/5.E-4_sub_tstep10_BC.h5','/t_2.0000/oz_ori');
u_re=h5read('../err_file/spline/5.E-4_sub_tstep10_BC.h5','/t_2.0000/u');
u_os=h5read('../AB2_result.MARCC/HIT_256^3_decay_5.E-4_AB2_dp_x0_16_nx0_32_sub.h5','/t_2.0000/u_sub');
v_re=h5read('../err_file/spline/5.E-4_sub_tstep10_BC.h5','/t_2.0000/v');
v_os=h5read('../AB2_result.MARCC/HIT_256^3_decay_5.E-4_AB2_dp_x0_16_nx0_32_sub.h5','/t_2.0000/v_sub');

dx=2*pi/256;
temp_re=zeros(33,33,34);
temp_os=zeros(33,33,34);
for i=1:33
    temp_re(i,:,:)=(v_re(i+1,:,:)-v_re(i,:,:))/dx;
    temp_os(i,:,:)=(v_os(i+1,:,:)-v_os(i,:,:))/dx;
end
for i=1:33
    temp_re(:,i,:)=temp_re(:,i,:)-(u_re(:,i+1,:)-u_re(:,i,:))/dx;
    temp_os(:,i,:)=temp_os(:,i,:)-(u_os(:,i+1,:)-u_os(:,i,:))/dx;
end

err=abs(u_re-u_os);
max(err(:))
idx=find(err==max(err(:)));
[d1,d2,d3]=ind2sub(size(err),idx);
d1=d1(1); d2=d2(1); d3=d3(1);

temp_re=squeeze(u_re(:,:,d3));
temp_os=squeeze(u_os(:,:,d3));
dx=2*pi/256;

figure;
[a,b]=contourf(temp_os(2:end-1,2:end-1),[-5:0.1:5],'LineStyle','none');% hold on;
temp=b.LevelList;
[a,b]=contourf(temp_os(2:end-1,2:end-1),temp,'LineStyle','none');
hold on;
contour(temp_re(2:end-1,2:end-1),temp,'k--')
caxis([-0.5,0.5])
clabel(a,b,[-4:1:4])
daspect([1 1 1])
xlabel('$i$','Interpreter','latex')
ylabel('$j$','Interpreter','latex')

font_size=10;
width=8; height=width*0.7;
set_figure(gca, gcf, width, height, font_size);
