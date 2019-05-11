clear all
clc

nx=32; ny=nx; nz=nx; t_start=1;
dt0=[4e-3,2e-3,1e-3,5e-4,2.5e-4];
sub_tstep=[2,3,4,5,6,7,8,9,10,20,50,100,1000];
t_span=[-4:5];
interp_scheme='spline';

temp=length(t_span);
u0=zeros(temp, nx+1, ny+2, nz+2);
v0=zeros(temp, nx+2, ny+1, nz+2);
w0=zeros(temp, nx+2, ny+2, nz+1);
p0=zeros(temp, nx+2, ny+2, nz+2);
u_star0=zeros(temp, nx+1, ny+2, nz+2);
v_star0=zeros(temp, nx+2, ny+1, nz+2);
w_star0=zeros(temp, nx+2, ny+2, nz+1);
dp0=zeros(temp, nx+2, ny+2, nz+2);

filename0(1)="../AB2_result.MARCC/HIT_256^3_decay_4.E-3_AB2_dp_x0_32_nx0_32_sub.h5";
filename0(2)="../AB2_result.MARCC/HIT_256^3_decay_2.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(3)="../AB2_result.MARCC/HIT_256^3_decay_1.E-3_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(4)="../AB2_result.MARCC/HIT_256^3_decay_5.E-4_AB2_dp_x0_16_nx0_32_sub.h5";
filename0(5)="../AB2_result.MARCC/HIT_256^3_decay_3.E-4_AB2_dp_x0_16_nx0_32_sub.h5";

for j=5:5
    idx=j; dt=dt0(idx);
    filename=char(filename0(idx));
    output_filename=[filename(1:end-3),'_',interp_scheme,filename(end-2:end)]
    t0=t_start+t_span*dt; 
    
    for i=1:length(t0)
        temp2=num2str(t0(i),'%.4f');
        if (t0(i)<1) temp2=temp2(2:end); end
        u0(i,:,:,:)=h5read(filename,['/t_',temp2,'/u_sub']);
        v0(i,:,:,:)=h5read(filename,['/t_',temp2,'/v_sub']);
        w0(i,:,:,:)=h5read(filename,['/t_',temp2,'/w_sub']);
        p0(i,:,:,:)=h5read(filename,['/t_',temp2,'/p_sub']);
        u_star0(i,:,:,:)=h5read(filename,['/t_',temp2,'/u_star_sub']);
        v_star0(i,:,:,:)=h5read(filename,['/t_',temp2,'/v_star_sub']);
        w_star0(i,:,:,:)=h5read(filename,['/t_',temp2,'/w_star_sub']);
        dp0(i,:,:,:)=h5read(filename,['/t_',temp2,'/dp_sub']);
    end
    
    for i=1:length(sub_tstep)
        t=linspace(t_start,t_start+dt,sub_tstep(i)+1); t=t(2:end);
        u=interp1(t0, u0, t, interp_scheme);
        v=interp1(t0, v0, t, interp_scheme);
        w=interp1(t0, w0, t, interp_scheme);
        p=interp1(t0, p0, t, interp_scheme);
        u_star=interp1(t0, u_star0, t, interp_scheme);
        v_star=interp1(t0, v_star0, t, interp_scheme);
        w_star=interp1(t0, w_star0, t, interp_scheme);
        dp=interp1(t0, dp0, t, interp_scheme);
        plot(t0, squeeze(dp0(:,5,6,20)),'-s', t, squeeze(dp(:,5,6,20)),'x');
        
        group_name=['sub_t_',num2str(sub_tstep(i))]
        h5create(output_filename, ['/',group_name,'/u_sub'], size(u));
        h5write(output_filename, ['/',group_name,'/u_sub'], u);
        h5create(output_filename, ['/',group_name,'/v_sub'], size(v));
        h5write(output_filename, ['/',group_name,'/v_sub'], v);
        h5create(output_filename, ['/',group_name,'/w_sub'], size(w));
        h5write(output_filename, ['/',group_name,'/w_sub'], w);
        h5create(output_filename, ['/',group_name,'/p_sub'], size(p));
        h5write(output_filename, ['/',group_name,'/p_sub'], p);
        h5create(output_filename, ['/',group_name,'/u_star_sub'], size(u_star));
        h5write(output_filename, ['/',group_name,'/u_star_sub'], u_star);
        h5create(output_filename, ['/',group_name,'/v_star_sub'], size(v_star));
        h5write(output_filename, ['/',group_name,'/v_star_sub'], v_star);
        h5create(output_filename, ['/',group_name,'/w_star_sub'], size(w_star));
        h5write(output_filename, ['/',group_name,'/w_star_sub'], w_star);
        h5create(output_filename, ['/',group_name,'/dp_sub'], size(dp));
        h5write(output_filename, ['/',group_name,'/dp_sub'], dp);
    end
end