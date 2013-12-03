function [alpha_x,alpha_y, fr]=get_response_AOD(data,f,xy_slopes,xy_k,cal,rate)

xy_k=-xy_k;
x=-data(1,:)/cal(5)*1e-6+1./(xy_slopes(1).*1e6).*data(4,:)./data(8,:);
%subplot(3,2,1)
%plot(-data(1,:)/cal(5)*1e-6)
%subplot(3,2,3)
%plot(1./(xy_slopes(1).*1e6).*data(4,:)./data(8,:))

y=-data(2,:)/cal(6)*1e-6+1./(xy_slopes(2).*1e6).*data(6,:)./data(8,:);
%subplot(3,2,2)
%plot(-data(2,:)/cal(6)*1e-6)
%subplot(3,2,4)
%plot(1./(xy_slopes(2).*1e6).*data(6,:)./data(8,:))

%subplot(3,2,5)
%plot(-data(1,:)/cal(5)*1e-6,1./(xy_slopes(1).*1e6).*data(4,:)./data(8,:))

%subplot(3,2,6)
%plot(-data(2,:)/cal(6)*1e-6,1./(xy_slopes(2).*1e6).*data(6,:)./data(8,:))

Fx=xy_k(1)./(xy_slopes(1)*1e6).*data(4,:)./data(8,:);
Fy=xy_k(2)./(xy_slopes(2)*1e6).*data(6,:)./data(8,:);
p=length(x);

alpha_x=fft(x)./fft(Fx);
alpha_y=fft(y)./fft(Fy);

alpha_x=alpha_x(1:p/2+1);
alpha_y=alpha_y(1:p/2+1);
fr=rate/p*([0:p/2]);
