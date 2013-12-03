function check_response
path='E:\Science\data\response_function\oocytes\2013-09-20\cellWT2_ves2_pbs1x\response_function_x=0.200_y=0.000_1000'
files=dir([path,filesep,'*deformation_response.mat']);
tic
pf=1.2;
p2=1
for j=1:length(files)
    load([path,filesep,files(j).name]);
    [alpha_x,alpha_y, fr]=get_response_AOD(squeeze(data),f,pf*xy_slope,p2*xy_k,cal,s_eff);
    %now pick the right value
    [a,b]=min(abs(f-fr));
    ax(j)=alpha_x(b);
    ay(j)=alpha_y(b);
    freq(j)=f;
end
G=1./(6*pi*1e-6*ax);
loglog(freq,abs(real(G)),freq,abs(imag(G)));
toc