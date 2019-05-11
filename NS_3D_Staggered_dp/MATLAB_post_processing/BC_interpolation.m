clear all
%clc

nx=32; ny=nx; nz=nx;
dt0=[4e-3,2e-3,1e-3,5e-4,2.5e-4];
idx=2; dt=dt0(idx);

t=[0.9:dt:2.1];
t_idx=find(t>=1 & t<=2);
step=5;
t0=t(1:step:end);
interp_scheme='spline';
snapshot_used=max(6,4);

temp=length(t);
u0=zeros(temp, nx+1, ny+2, nz+2);
v0=zeros(temp, nx+2, ny+1, nz+2);
w0=zeros(temp, nx+2, ny+2, nz+1);
p0=zeros(temp, nx+2, ny+2, nz+2);
u_star0=zeros(temp, nx+1, ny+2, nz+2);
v_star0=zeros(temp, nx+2, ny+1, nz+2);
w_star0=zeros(temp, nx+2, ny+2, nz+1);
dp0=zeros(temp, nx+2, ny+2, nz+2);

filename0(1)="../AB2_result.MARCC/HIT_256^3_decay_4.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(2)="../AB2_result.MARCC/HIT_256^3_decay_2.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(3)="../AB2_result.MARCC/HIT_256^3_decay_1.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(4)="../AB2_result.MARCC/HIT_256^3_decay_5.E-4_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(5)="../AB2_result.MARCC/HIT_256^3_decay_3.E-4_AB2_dp_x0_16_nx0_32_sub.h5";

filename=char(filename0(idx));
output_filename=[filename(1:end-3),'_BC_',interp_scheme,filename(end-2:end)]

for i=1:length(t)
    temp2=num2str(t(i),'%.4f');
    if (t(i)<1) temp2=temp2(2:end); end
%     u0(i,:,:,:)=h5read(filename,['/t_',temp2,'/u_sub']);
%     v0(i,:,:,:)=h5read(filename,['/t_',temp2,'/v_sub']);
%     w0(i,:,:,:)=h5read(filename,['/t_',temp2,'/w_sub']);
%     p0(i,:,:,:)=h5read(filename,['/t_',temp2,'/p_sub']);
    u_star0(i,:,:,:)=h5read(filename,['/t_',temp2,'/u_star_sub']);
    v_star0(i,:,:,:)=h5read(filename,['/t_',temp2,'/v_star_sub']);
    w_star0(i,:,:,:)=h5read(filename,['/t_',temp2,'/w_star_sub']);
    dp0(i,:,:,:)=h5read(filename,['/t_',temp2,'/dp_sub']);
end

u_star_base=u_star0(1:step:end,:,:,:); v_star_base=v_star0(1:step:end,:,:,:);
w_star_base=w_star0(1:step:end,:,:,:); dp_base=dp0(1:step:end,:,:,:);

u_star=zeros(length(t_idx), nx+1, ny+2, nz+2); v_star=zeros(length(t_idx), nx+2, ny+1, nz+2);
w_star=zeros(length(t_idx), nx+2, ny+2, nz+1); dp=zeros(length(t_idx), nx+2, ny+2, nz+2);

for i=1:length(t_idx)
    t_int=t(t_idx(i));
%     if (t_int<=1.9)
    temp=find(t0<=t_int); temp=temp(end-(snapshot_used/2)+1); temp=temp:temp+snapshot_used-1;
    % u=interp1(t0, u0(1:step:end,:,:,:), t(t_idx), interp_scheme);
    % v=interp1(t0, v0, t, interp_scheme);
    % w=interp1(t0, w0, t, interp_scheme);
    % p=interp1(t0, p0, t, interp_scheme);
    u_star(i,:,:,:)=interp1(t0(temp), u_star_base(temp,:,:,:), t_int, interp_scheme);
    v_star(i,:,:,:)=interp1(t0(temp), v_star_base(temp,:,:,:), t_int, interp_scheme);
    w_star(i,:,:,:)=interp1(t0(temp), w_star_base(temp,:,:,:), t_int, interp_scheme);
    dp(i,:,:,:)=interp1(t0(temp), dp_base(temp,:,:,:), t_int, interp_scheme);
%     else
%         u_star(i,:,:,:)=u_star0(t_idx(i),:,:,:);
%         v_star(i,:,:,:)=v_star0(t_idx(i),:,:,:);
%         w_star(i,:,:,:)=w_star0(t_idx(i),:,:,:);
%         dp(i,:,:,:)=dp0(t_idx(i),:,:,:);
%     end
    
%     temp2=num2str(t_int,'%.4f');
%     if (t_int<1) temp2=temp2(2:end); end
%     group_name=['t_',num2str(temp2)]
% %     h5create(output_filename, ['/',group_name,'/u_sub'], size(u));
% %     h5write(output_filename, ['/',group_name,'/u_sub'], u);
% %     h5create(output_filename, ['/',group_name,'/v_sub'], size(v));
% %     h5write(output_filename, ['/',group_name,'/v_sub'], v);
% %     h5create(output_filename, ['/',group_name,'/w_sub'], size(w));
% %     h5write(output_filename, ['/',group_name,'/w_sub'], w);
% %     h5create(output_filename, ['/',group_name,'/p_sub'], size(p));
% %     h5write(output_filename, ['/',group_name,'/p_sub'], p);
%     temp1=squeeze(u_star(i,:,:,:));
%     temp2=squeeze(v_star(i,:,:,:));
%     temp3=squeeze(w_star(i,:,:,:));
%     temp4=squeeze(dp(i,:,:,:));
%     h5create(output_filename, ['/',group_name,'/u_star_sub'], size(temp1));
%     h5write(output_filename, ['/',group_name,'/u_star_sub'], temp1);
%     h5create(output_filename, ['/',group_name,'/v_star_sub'], size(temp2));
%     h5write(output_filename, ['/',group_name,'/v_star_sub'], temp2);
%     h5create(output_filename, ['/',group_name,'/w_star_sub'], size(temp3));
%     h5write(output_filename, ['/',group_name,'/w_star_sub'], temp3);
%     h5create(output_filename, ['/',group_name,'/dp_sub'], size(temp4));
%     h5write(output_filename, ['/',group_name,'/dp_sub'], temp4);
end

plot(t0, squeeze(dp0(1:step:end,5,6,20)),'-s', t(t_idx), squeeze(dp(:,5,6,20)),'x');
%%
t_start=1; t_end=2;
color=[0,0.45,0.74;0.85,0.33,0.10;0.93,0.69,0.13];
dt=[4e-3,2e-3,1e-3,5e-4,2.5e-4]; 
t1=[t_start:dt(1):t_end]; t2=[t_start:dt(2):t_end]; t3=[t_start:dt(3):t_end]; t4=[t_start:dt(4):t_end]; t5=[t_start:dt(5):t_end];

u_star(:,2:end-1,3:end-2,3:end-2)=nan;
temp0=u_star0(t_idx,:,:,:)-u_star; temp0(:,2:end-1,3:end-2,3:end-2)=nan; temp1=max(max(max(abs(temp0),[],4),[],3),[],2)/rms(u_star(~isnan(u_star)));
v_star(:,3:end-2,2:end-1,3:end-2)=nan;
temp0=v_star0(t_idx,:,:,:)-v_star; temp0(:,3:end-2,2:end-1,3:end-2)=nan; temp2=max(max(max(abs(temp0),[],4),[],3),[],2)/rms(v_star(~isnan(v_star)));
w_star(:,3:end-2,3:end-2,2:end-1)=nan;
temp0=w_star0(t_idx,:,:,:)-w_star; temp0(:,3:end-2,3:end-2,2:end-1)=nan; temp3=max(max(max(abs(temp0),[],4),[],3),[],2)/rms(w_star(~isnan(w_star)));
dp(:,3:end-2,3:end-2,3:end-2)=nan;
temp0=dp0(t_idx,:,:,:)-dp; temp0(:,3:end-2,3:end-2,3:end-2)=nan; temp4=max(max(max(abs(temp0),[],4),[],3),[],2)/rms(dp(~isnan(dp)));

[max(temp1) max(temp2) max(temp3) max(temp4)]

figure
plot(t2,temp1,'-')
hold on;
plot(t2,temp2,'--')
plot(t2,temp3,'-.')
plot(t2,temp4,'-k')
legend('u','v','w','p')
xlabel('$t$','Interpreter','latex')
ylabel('$\epsilon_{\varphi,L_\infty}$','Interpreter','latex')

