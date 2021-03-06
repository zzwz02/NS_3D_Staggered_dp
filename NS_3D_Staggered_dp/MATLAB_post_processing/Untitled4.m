clear all
%close all
%clc

t_start=1; t_end=2;
dt=[4e-3, 2e-3, 1e-3, 5e-4];
dt_string={'4.E-3','2.E-3','1.E-3','5.E-4'};

err_u_star_max=zeros(length(length(dt)),1);
err_RHS_max=err_u_star_max; err_dp_max=err_u_star_max; err_p_max=err_u_star_max; err_u_max=err_u_star_max;

for T=1:20
    for idx=1:length(dt)
        t=[t_start:dt(idx):t_end];
        
        resim_file=['../err_file/err_vel_AB2_',dt_string{idx},'_p_N_Ustar_TOffset_restart_noise1.E-4.h5'];
        bigDNS_file=['../',h5readatt(resim_file,'/','big_DNS_file')];
        slice=16;
        
        temp=h5info(resim_file);
        %if (idx==1)
            gName=temp.Groups(T+1).Name;
        %end
        u_star=h5read(resim_file,[gName,'/u_star_sub']); %u=u(:,2:end-1,2:end-1);
        RHS=h5read(resim_file,[gName,'/RHS_poisson_sub']); RHS=RHS(2:end-1,2:end-1,2:end-1);
        dp=h5read(resim_file,[gName,'/dp_sub']); %dp=dp-max(dp(:));
        p=h5read(resim_file,[gName,'/p_sub']); %p=p-max(p(:));
        u=h5read(resim_file,[gName,'/u_sub']); %u=u(:,2:end-1,2:end-1);
        
        u_star_sub=h5read(bigDNS_file,[gName,'/u_star_sub']); %u_sub=u_sub(:,2:end-1,2:end-1);
        RHS_sub=h5read(bigDNS_file,[gName,'/RHS_poisson_sub']); RHS_sub=RHS_sub(2:end-1,2:end-1,2:end-1);
        dp_sub=h5read(bigDNS_file,[gName,'/dp_sub']); %dp_sub=dp_sub-max(dp_sub(:));
        p_sub=h5read(bigDNS_file,[gName,'/p_sub']); %p_sub=p_sub-max(p_sub(:));
        u_sub=h5read(bigDNS_file,[gName,'/u_sub']); %u=u(:,2:end-1,2:end-1);
        
        err_u_star=(u_star-u_star_sub)/rms(u_star_sub(:));
        err_RHS=(RHS-RHS_sub)/rms(RHS_sub(:));
        err_dp=(dp-dp_sub)/rms(dp_sub(:));
        err_p=(p-p_sub)/rms(p_sub(:));
        err_u=(u-u_sub)/rms(u_sub(:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp1=u_star(2:end-1,2:end-1,2:end-1); temp2=u_star_sub(2:end-1,2:end-1,2:end-1);
        err_u_star_in=(temp1-temp2)/rms(temp2(:));
        
        temp1=RHS(2:end-1,2:end-1,2:end-1); temp2=RHS_sub(2:end-1,2:end-1,2:end-1);
        err_RHS_in=(temp1-temp2)/rms(temp2(:));
        
        temp1=dp(2:end-1,2:end-1,2:end-1); temp2=dp_sub(2:end-1,2:end-1,2:end-1);
        err_dp_in=(temp1-temp2)/rms(temp2(:));
        
        temp1=p(2:end-1,2:end-1,2:end-1); temp2=p_sub(2:end-1,2:end-1,2:end-1);
        err_p_in=(temp1-temp2)/rms(temp2(:));
        
        temp1=u(2:end-1,2:end-1,2:end-1); temp2=u_sub(2:end-1,2:end-1,2:end-1);
        err_u_in=(temp1-temp2)/rms(temp2(:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp1=u_star; temp2=u_star_sub; temp1(2:end-1,2:end-1,2:end-1)=nan; temp2(2:end-1,2:end-1,2:end-1)=nan;
        temp1=temp1(~isnan(temp1(:)));temp2=temp2(~isnan(temp2(:)));
        err_u_star_bc=(temp1-temp2)/rms(temp2(:));
        
        temp1=RHS; temp2=RHS_sub; temp1(2:end-1,2:end-1,2:end-1)=nan; temp2(2:end-1,2:end-1,2:end-1)=nan;
        temp1=temp1(~isnan(temp1(:)));temp2=temp2(~isnan(temp2(:)));
        err_RHS_bc=(temp1-temp2)/rms(temp2(:));
        
        temp1=dp; temp2=dp_sub; temp1(2:end-1,2:end-1,2:end-1)=nan; temp2(2:end-1,2:end-1,2:end-1)=nan;
        temp1=temp1(~isnan(temp1(:)));temp2=temp2(~isnan(temp2(:)));
        err_dp_bc=(temp1-temp2)/rms(temp2(:));
        
        temp1=p; temp2=p_sub; temp1(2:end-1,2:end-1,2:end-1)=nan; temp2(2:end-1,2:end-1,2:end-1)=nan;
        temp1=temp1(~isnan(temp1(:)));temp2=temp2(~isnan(temp2(:)));
        err_p_bc=(temp1-temp2)/rms(temp2(:));
        
        temp1=u; temp2=u_sub; temp1(2:end-1,2:end-1,2:end-1)=nan; temp2(2:end-1,2:end-1,2:end-1)=nan;
        temp1=temp1(~isnan(temp1(:)));temp2=temp2(~isnan(temp2(:)));
        err_u_bc=(temp1-temp2)/rms(temp2(:));
        
        err_u_star_max(idx)=max(abs(err_u_star(:)));
        err_RHS_max(idx)=max(abs(err_RHS(:)));
        err_dp_max(idx)=max(abs(err_dp(:)));
        err_p_max(idx)=max(abs(err_p(:)));
        err_u_max(idx)=max(abs(err_u(:)));
        
        err_u_star_in_max(idx)=max(abs(err_u_star_in(:)));
        err_RHS_in_max(idx)=max(abs(err_RHS_in(:)));
        err_dp_in_max(idx)=max(abs(err_dp_in(:)));
        err_p_in_max(idx)=max(abs(err_p_in(:)));
        err_u_in_max(idx)=max(abs(err_u_in(:)));
        
        err_u_star_bc_max(idx)=max(abs(err_u_star_bc(:)));
        err_RHS_bc_max(idx)=max(abs(err_RHS_bc(:)));
        err_dp_bc_max(idx)=max(abs(err_dp_bc(:)));
        err_p_bc_max(idx)=max(abs(err_p_bc(:)));
        err_u_bc_max(idx)=max(abs(err_u_bc(:)));
    end
    
    clear temp1 temp2
    temp1(1,:)=polyfit(log(dt),log(err_u_star_max),1);
    temp1(2,:)=polyfit(log(dt),log(err_RHS_max),1);
    temp1(3,:)=polyfit(log(dt),log(err_dp_max),1);
    temp1(4,:)=polyfit(log(dt),log(err_p_max),1);
    temp1(5,:)=polyfit(log(dt),log(err_u_max),1);
    
    temp2(1,:)=polyfit(log(dt),log(err_u_star_in_max),1);
    temp2(2,:)=polyfit(log(dt),log(err_RHS_in_max),1);
    temp2(3,:)=polyfit(log(dt),log(err_dp_in_max),1);
    temp2(4,:)=polyfit(log(dt),log(err_p_in_max),1);
    temp2(5,:)=polyfit(log(dt),log(err_u_in_max),1);
    
    temp3(1,:)=polyfit(log(dt),log(err_u_star_bc_max),1);
    temp3(2,:)=polyfit(log(dt),log(err_RHS_bc_max),1);
    temp3(3,:)=polyfit(log(dt),log(err_dp_bc_max),1);
    temp3(4,:)=polyfit(log(dt),log(err_p_bc_max),1);
    temp3(5,:)=polyfit(log(dt),log(err_u_bc_max),1);
    disp([num2str(T),'     u*        RHS       dp        p        u'])
    disp([temp1(:,1)'; temp2(:,1)'; temp3(:,1)'])
    
    subplot(4,5,T)
    loglog(dt,err_u_star_max,'k^',dt,err_RHS_max,'k+',dt,err_dp_max,'kd',dt,err_p_max,'ks',dt,err_u_max,'kv');
    %legend('u*','RHS','dp','p','u')
    hold on
    plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^2*0.7e2,'-','Color',[0.47,0.67,0.19])
    plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^1*3e3,'--','Color',[0.47,0.67,0.19])
    plot([dt(1)*1.1 dt(end)*0.9],[dt(1)*1.1 dt(end)*0.9].^1*5e0,'--','Color',[0.47,0.67,0.19])
    title(num2str(temp1(:,1)','%.2f '))
end