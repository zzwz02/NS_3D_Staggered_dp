clear all
clc

nx=32; ny=nx; nz=nx; t_start=1;
u0=zeros(nx+1, ny+2, nz+2, 4);
v0=zeros(nx+2, ny+1, nz+2, 4);
w0=zeros(nx+2, ny+2, nz+1, 4);
dp0=zeros(nx+2, ny+2, nz+2, 4);
filename0(1)="../AB2_result.MARCC/HIT_256^3_decay_4.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(2)="../AB2_result.MARCC/HIT_256^3_decay_2.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(3)="../AB2_result.MARCC/HIT_256^3_decay_1.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(4)="../AB2_result.MARCC/HIT_256^3_decay_5.E-4_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(5)="../AB2_result.MARCC/HIT_256^3_decay_3.E-4_AB2_dp_x0_16_nx0_32_sub.h5";
dt0=[4e-3,2e-3,1e-3,5e-4,2.5e-4];

for j=1:length(dt0)
    idx=j; dt=dt0(idx); filename=char(filename0(idx)); output_filename=[filename(1:end-3),'_pchip',filename(end-2:end)]
    t0=t_start+[-1:2]*dt; sub_tstep=[2,5,10,20,50,100];
    
    for i=1:4
        temp2=num2str(t0(i),'%.4f');
        if (t0(i)<1) temp2=temp2(2:end); end
        u0(:,:,:,i)=h5read(filename,['/t_',temp2,'/u_sub']);
        v0(:,:,:,i)=h5read(filename,['/t_',temp2,'/v_sub']);
        w0(:,:,:,i)=h5read(filename,['/t_',temp2,'/w_sub']);
        dp0(:,:,:,i)=h5read(filename,['/t_',temp2,'/dp_sub']);
    end
    
    for i=1:length(sub_tstep)
        t=linspace(t_start,t_start+dt,sub_tstep(i)+1); t=t(2:end);
        u=pchip(t0, u0, t);
        v=pchip(t0, v0, t);
        w=pchip(t0, w0, t);
        dp=pchip(t0, dp0, t);
        %plot(t0, squeeze(dp0(5,6,20,:)),'-s', t, squeeze(dp(5,6,20,:)),'x');
        
        group_name=['sub_t_',num2str(sub_tstep(i))]
        h5create(output_filename, ['/',group_name,'/u'], size(u));
        h5write(output_filename, ['/',group_name,'/u'], u);
        h5create(output_filename, ['/',group_name,'/v'], size(v));
        h5write(output_filename, ['/',group_name,'/v'], v);
        h5create(output_filename, ['/',group_name,'/w'], size(w));
        h5write(output_filename, ['/',group_name,'/w'], w);
        h5create(output_filename, ['/',group_name,'/dp'], size(dp));
        h5write(output_filename, ['/',group_name,'/dp'], dp);
    end
end