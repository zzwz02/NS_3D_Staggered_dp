clear all
close all
clc

color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13;0.49,0.18,0.56];

dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; sub_tstep=[1,2,3,4,5,6,7,8,9,10,20,50,100];
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2]; t5=[1:dt(5):2];

err{1} = import_err_file('../err_file/spline/base_err_vel_AB2_2.E-3_N_Ustar_TOffset.txt');
for i=2:length(sub_tstep)
    err{i} = import_err_file(['../err_file/spline/base_err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep',num2str(sub_tstep(i)),'.txt']);
end 
err{length(sub_tstep)+1} = import_err_file('../err_file/spline/base_err_vel_AB2_2.E-3_N_Ustar_TOffset_restart.txt');
%%

tt=t2;
figure(1)
semilogy(tt,err{1}.rel_err_p)
hold on;
plot(tt,err{2}.rel_err_p)
plot(tt,err{5}.rel_err_p)
plot(tt,err{10}.rel_err_p)
plot(tt,err{end}.rel_err_p,'--')

err_u_max=zeros(length(err),1); err_v_max=err_u_max; err_w_max=err_u_max; err_p_max=err_u_max; 
crit=1.02;
for i=1:length(err)
    err_u_max(i)=max(err{i}.rel_err_u(tt>1.02));
    err_v_max(i)=max(err{i}.rel_err_v(tt>1.02));
    err_w_max(i)=max(err{i}.rel_err_w(tt>1.02));
    err_p_max(i)=max(err{i}.rel_err_p(tt>1.02));
end

figure(2)
loglog(sub_tstep,err_u_max(1:end-1),'^','Color',color(1,:))
hold on;
loglog(sub_tstep,err_v_max(1:end-1),'s','Color',color(2,:))
loglog(sub_tstep,err_w_max(1:end-1),'o','Color',color(3,:))
loglog(sub_tstep,err_p_max(1:end-1),'+','Color',color(4,:))

loglog([sub_tstep(1) sub_tstep(end)],[err_u_max(end) err_u_max(end)],'-','Color',color(1,:))
loglog([sub_tstep(1) sub_tstep(end)],[err_v_max(end) err_v_max(end)],'--','Color',color(2,:))
loglog([sub_tstep(1) sub_tstep(end)],[err_w_max(end) err_w_max(end)],'-.','Color',color(3,:))
loglog([sub_tstep(1) sub_tstep(end)],[err_p_max(end) err_p_max(end)],':','Color',color(4,:))

loglog([sub_tstep(1)*0.9 sub_tstep(end)*1.1],[sub_tstep(1)*0.9 sub_tstep(end)*1.1].^-2*3e-3)

legend('u','v','w','p')
xlabel('number of sub-time step $l$','Interpreter','latex'); ylabel('$\epsilon_{\varphi,L_\infty}$ (compared with reference solut.)','Interpreter','latex')

temp(1,:)=abs(err_u_max(1:end-1)-err_u_max(end));
temp(2,:)=abs(err_v_max(1:end-1)-err_v_max(end));
temp(3,:)=abs(err_w_max(1:end-1)-err_w_max(end));
temp(4,:)=abs(err_p_max(1:end-1)-err_p_max(end));

figure(3)
loglog(sub_tstep,temp(1,:),'-+',sub_tstep,temp(2,:),'--+',sub_tstep,temp(3,:),'-.+',sub_tstep,temp(4,:),':+')

%%
clear all

color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13;0.49,0.18,0.56];

dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; sub_tstep=[1,2,3,4,5,6,7,8,9,10,20,50,100];
t1=[1:dt(1):2]; t2=[1:dt(2):2]; t3=[1:dt(3):2]; t4=[1:dt(4):2]; t5=[1:dt(5):2];

err{1} = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset.txt');
for i=2:length(sub_tstep)
    err{i} = import_err_file(['../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_sub_tstep',num2str(sub_tstep(i)),'.txt']);
end 
%err{length(sub_tstep)+1} = import_err_file('../err_file/spline/err_vel_AB2_2.E-3_N_Ustar_TOffset_restart.txt');

tt=t2;
pcolor=color(1,:);
figure(4)
semilogy(tt,err{1}.rel_err_p,'--')
hold on;
semilogy(tt,err{2}.rel_err_p)
semilogy(tt,err{5}.rel_err_p)
semilogy(tt,err{10}.rel_err_p)
semilogy(tt,err{11}.rel_err_p)

err_u_max=zeros(length(err),1); err_v_max=err_u_max; err_w_max=err_u_max; err_p_max=err_u_max; 
crit=1.02;
for i=1:length(err)
    err_u_max(i)=max(err{i}.rel_err_u(tt>crit));
    err_v_max(i)=max(err{i}.rel_err_v(tt>crit));
    err_w_max(i)=max(err{i}.rel_err_w(tt>crit));
    err_p_max(i)=max(err{i}.rel_err_p(tt>crit));
end

figure(5)
loglog(sub_tstep,err_u_max,'^','Color',color(1,:))
hold on;
loglog(sub_tstep,err_v_max,'s','Color',color(2,:))
loglog(sub_tstep,err_w_max,'o','Color',color(3,:))
loglog(sub_tstep,err_p_max,'+','Color',color(4,:))

loglog([1*0.9 10*1.1],[1*0.9 10*1.1].^-2*4.5e-4,'k--')
xlim([0.9, 103])
polyfit(log(sub_tstep(1:4)),log(err_u_max(1:4))',1)

legend('u','v','w','p','slope=2')
xlabel('number of sub-time step $l$','Interpreter','latex'); ylabel('$\epsilon_{\varphi,L_\infty}$ (compared with AB2 solut.)','Interpreter','latex')
