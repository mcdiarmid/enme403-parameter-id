clc
beep off
clear all
close all
format short g
set(0,'units','centimeters','defaultfigureposition',[60 60 420,330])
figure(1201)

% back to the Creatinine model
t=0:600;
TT=0:60:600;
UX=zeros(size(t));UX(62)=4;
V=40;
k=0.004;
UN=0.01;
C=exp(-k*t).*(UN/V/k+cumtrapz(t,exp(k*t).*(UX+UN)/V));
CC=C(TT+1);
CC=[0.0617362286206992,0.0578591828680488,0.138958525405719,0.126835540580120,0.121588113934300,0.105473590752212,0.0939101320790592,0.0905332846505024,0.0791227372065423,0.0721591272672233,0.0715330339740876;];
CC0=mean(CC(1:2));

I=0;

%determine the objective surface!
for k=0:0.0001:0.008
    I=I+1;J=0;
    for UN=0:0.0002:0.02
        J=J+1;
        CCt=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(UX+UN)/V));
        psi(I,J)=norm(CCt(TT+1)-CC);
    end
end

% contour(0:0.0001:0.008,0:0.0002:0.02,log(psi'),50);hold all

alpha=0.1;

KUN=[0.00, 0.0200]; %I.C.
dKUN=[0.000004 0.00001]; %perturbations

for III=1:1000
    CCt=exp(-KUN(III,1)*t).*(CC0+cumtrapz(t,exp(KUN(III,1)*t).*(UX+KUN(III,2))/V)); %model simulation
    psi2(III)=norm(CCt(TT+1)-CC);  % The 2-norm of the model simulation (CCt) at the sample times (TT+1) minus the raw data
    
    CCtk=exp(-(KUN(III,1)+dKUN(1))*t).*(CC0+cumtrapz(t,exp((KUN(III,1)+dKUN(1))*t).*(UX+KUN(III,2))/V));
    nCCtk=norm(CCtk(TT+1)-CC);
    
    CCtu=exp(-KUN(III,1)*t).*(CC0+cumtrapz(t,exp(KUN(III,1)*t).*(UX+(KUN(III,2)+dKUN(2)))/V));
    nCCtu=norm(CCtu(TT+1)-CC);
    
    % the Jacobian
    JJ=[(nCCtk-psi2(III))/dKUN(1) (nCCtu-psi2(III))/dKUN(2)];
    
    % the gradient descent bit
    KUN(III+1,:)=KUN(III,:)-alpha*JJ;
end

plot(KUN(:,1),KUN(:,2),'.-k');hold all
plot(KUN(1,1),KUN(1,2),'ok')
plot(0.00207, 0.002271,'+k')

xlabel('\it\bfk')
ylabel('\it\bfU_N')

figure(1202)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,0:0.0001:0.008,0:0.0002:0.02,psi',50);hold all
alpha=0.000001;         % make more stable

KUN=[0.000, 0.02];
dKUN=[0.000004 0.00001];

for III=1:20000
    CCt=exp(-KUN(III,1)*t).*(CC0+cumtrapz(t,exp(KUN(III,1)*t).*(UX+KUN(III,2))/V));
    psi2(III)=norm(CCt(TT+1)-CC);
    
    CCtk=exp(-(KUN(III,1)+dKUN(1))*t).*(CC0+cumtrapz(t,exp((KUN(III,1)+dKUN(1))*t).*(UX+KUN(III,2))/V));
    nCCtk=norm(CCtk(TT+1)-CC);
    
    CCtu=exp(-KUN(III,1)*t).*(CC0+cumtrapz(t,exp(KUN(III,1)*t).*(UX+(KUN(III,2)+dKUN(2)))/V));
    nCCtu=norm(CCtu(TT+1)-CC);
    
    JJ=[(nCCtk-psi2(III))/dKUN(1) (nCCtu-psi2(III))/dKUN(2)];
    KUN(III+1,:)=KUN(III,:)-alpha*JJ;
end

plot(h1,KUN(:,1),KUN(:,2),'.-k');hold all
plot(h1,KUN(1,1),KUN(1,2),'ok')
plot(h1,0.004, 0.01,'+r','markersize',10,'linewidth',2)

xlabel(h1,'\it\bfk')
ylabel(h1,'\it\bfU_N')
plot(h2,psi2,'.-k'); 
set(h2,'Yscale','log','ylim',[.005 1])
ylabel('\it\bf\psi')
xlabel('\bfiterations')

figure(1203)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,0:0.0001:0.008,0:0.0002:0.02,psi',50);hold all
alpha=0.000001;

KUN=[0.008, 0.015];
dKUN=[0.000004 0.00001];

for III=1:20000
    CCt=exp(-KUN(III,1)*t).*(CC0+cumtrapz(t,exp(KUN(III,1)*t).*(UX+KUN(III,2))/V));
    psi2(III)=norm(CCt(TT+1)-CC);
    
    CCtk=exp(-(KUN(III,1)+dKUN(1))*t).*(CC0+cumtrapz(t,exp((KUN(III,1)+dKUN(1))*t).*(UX+KUN(III,2))/V));
    nCCtk=norm(CCtk(TT+1)-CC);
    
    CCtu=exp(-KUN(III,1)*t).*(CC0+cumtrapz(t,exp(KUN(III,1)*t).*(UX+(KUN(III,2)+dKUN(2)))/V));
    nCCtu=norm(CCtu(TT+1)-CC);
    
    JJ=[(nCCtk-psi2(III))/dKUN(1) (nCCtu-psi2(III))/dKUN(2)];
    
    KUN(III+1,:)=KUN(III,:)-alpha*JJ;
    if III==12000
        figure(1204)
        plot(TT,CC,'+k',t,CCt,'b');hold all
    end
end
plot(t,CCt,'g');hold all
xlabel('bfTime [min]')
ylabel('\bfCreatinine [mmol.L^-^1]')
legend('data','optima','final');

figure(1203)
plot(h1,KUN(:,1),KUN(:,2),'.-k');hold all
plot(h1,KUN(1,1),KUN(1,2),'ok')
plot(h1,0.004, 0.01,'+r','markersize',10,'linewidth',2)

xlabel(h1,'\it\bfk')
ylabel(h1,'\it\bfU_N')
plot(h2,psi2,'.-k'); 
set(h2,'Yscale','log','ylim',[.005 1])
ylabel('\it\bf\psi')
xlabel('\bfiterations')
