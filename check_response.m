function check_response
path='E:\Science\data\response_function\oocytes\2013-09-20\cellWT10_ves1_pbs1x_slide2\response_function_x=0.200_y=0.000_1000'
files=dir([path,filesep,'*deformation_response.mat']);
tic
pf=1;
p2=1
for j=1:length(files)
    load([path,filesep,files(j).name]);
    [alpha_x,alpha_y, fr]=get_response_AOD(squeeze(data),f,pf*xy_slope,p2*xy_k,cal,s_eff);
    %now pick the right value
    [a,b]=min(abs(f-fr));
    ax(j)=alpha_x(b);
    ay(j)=alpha_y(b);
    freq(j)=f;
    data=squeeze(data);
    %and here I use the direct fit with the cosine and sine
    x=-data(1,:)/cal(5)*1e-6+1./(xy_slope(1).*1e6).*data(4,:)./data(8,:);
    Fx=xy_k(1)./(xy_slope(1)*1e6).*data(4,:)./data(8,:);
    temp=[0:length(x)-1]/s_eff;
    f_sin_f=fittype('a*sin(2*pi*f*t+phase)','independent','t','problem','f');
    [ff,goff]=fit(temp',Fx'-mean(Fx),f_sin_f,'problem',f);
    phase=ff.phase;


    f_sin_p=fittype('a*sin(2*pi*f*t+phase)+b*cos(2*pi*f*t+phase)','independent','t','problem',{'f','phase'});
    [fp,gofp]=fit(temp',x'-mean(x),f_sin_p,'problem',{f,phase});
    clear i;
    response_recal(j)=fp.a/ff.a+i*(fp.b/ff.a);
end
G=1./(6*pi*1e-6*ax);
G2=1./(6*pi*1e-6*response_recal);
loglog(freq,abs(real(G)),freq,abs(imag(G)));
hold on
loglog(freq,abs(imag(G)),'g')
loglog(freq,abs(real(G2)),'ro')
loglog(freq,abs(imag(G2)),'rx')
hold off

toc