clc
beep off
clear all
close all
format short g

% set(0,'units','centimeters','defaultfigureposition',[60 60 450 330])

t=(0:0.001:15)';
tt=(0:0.5:15)';
phi1=[ones(10001,1);zeros(5000,1)];
phi2=[zeros(10001,1);ones(5000,1)];


%       This commented section generates parent model and data - this woul
%       dnt happen in treal life....
a1=-0.001;
a2=-0.1;
ud=zeros(size(t));

for i=1:100
    ud=0+max(0,min(100, cumtrapz(t, 10*exp(-0.02*ud).*phi1+a1*phi1.*ud.^2+a2*phi2.*ud.^2)));
end

% ud_data=ud(1:50:1501)+randn(size(tt));
% ud_data(1)=0;

 ud_data=[0;6.14205297110053;8.14698995214801;13.7935478959043;16.8747933990616;19.7278135528891;22.9850537558254;25.7066805523054;28.4221671848922;30.8586498167937;32.0985010997414;33.6995084341626;37.5244501621723;38.0057740141821;41.1714163692344;40.2393841481550;41.8614105856496;45.0632483894187;46.1388684642106;44.7021780983235;46.8553968888880;15.4380634251729;9.60126348997106;4.56506174740374;6.14444844788981;4.35246346210841;3.34316635514530;2.47025569416698;1.40962409486031;2.21841570029254;1.66286180315103];
figure(811)
set(gcf,'position',[60 60 450 330])
plot(tt,ud_data,'+r');
xlabel('\bfTime [s]')
ylabel('\bfSpeed [m/s]')


figure(812)
set(gcf,'position',[60 60 680 330])
udG=interp1(tt,ud_data,t);

for i=1:50

    ia1=cumtrapz(t,phi1.*udG.*udG);
    ia2=cumtrapz(t,phi2.*udG.*udG);
    iFm=cumtrapz(t,phi1.*10.*exp(-0.02*udG));

    A=[ia1(1:500:15001) ia2(1:500:15001)];
    b=ud_data-ud_data(1)-iFm(1:500:15001);

    aa=A\b;
    
    for j=1:100
        udG=0+max(0,min(100, cumtrapz(t, 10*exp(-0.02*udG).*phi1+aa(1)*phi1.*udG.^2+aa(2)*phi2.*udG.^2)));
    end
    psis(i)=norm(ud_data-udG(1:500:15001));
    
   subplot(1,3,1:2);plot(t,udG,'b',tt,ud_data,'+r');xlabel('\bfTime [s]');ylabel('\bfSpeed [m/s]');drawnow;
   subplot(1,3,3);plot(tt,ud_data-udG(1:500:15001),'+r',tt,sort(ud_data-udG(1:500:15001)),'-b');xlabel('\bfTime [s]');title('\bfResiduals');axis([0 15 -3 3]);
   pause(0.0);
    
    
end

figure
plot(psis)
ylabel('\bf\it\psi')
xlabel('\bfiterations')


















