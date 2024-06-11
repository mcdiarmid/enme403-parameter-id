clc
beep off
clear all
close all
format short g
set(0,'units','centimeters','defaultfigureposition',[60 60 420,330])

t=0:600;
TT=0:60:600;
UX=zeros(size(t));UX(62)=4;
V=40;
k=0.004;
UN=0.01;
C=exp(-k*t).*(UN/V/k+cumtrapz(t,exp(k*t).*(UX+UN)/V));
CC=C(TT+1);
CC=[0.0617362286206992,0.0578591828680488,0.138958525405719,0.126835540580120,0.121588113934300,0.105473590752212,0.0939101320790592,0.0905332846505024,0.0791227372065423,0.0721591272672233,0.0715330339740876;];
% CC(1:2)=mean(CC(1:2));
CC0=mean(CC(1:2));

I=0;

%determine the objective surface!
for k=-.002:0.0001:0.008
    I=I+1;J=0;
    for UN=-.01:0.0002:0.02
        J=J+1;
        CCt=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(UX+UN)/V));
        psi(I,J)=norm(CCt(TT+1)-CC);
    end
end
figure(1301)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

alpha=1;
KUN=[0.00; 0.0200]; %I.C.
dKUN=[0.000004 0.00001]; %perturbations

for III=1:100
    CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
    psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
    psis(III)=norm(psi2(III,:));
    
    CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
    nCCtk=CCtk(TT+1)-CC;
    
    CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
    nCCtu=CCtu(TT+1)-CC;
    
    JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
    
    % this time with added Gauss Newton goodness:
    KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ)^-1*JJ'*psi2(III,1:11)';
end

plot(h1,KUN(1,:),KUN(2,:),'.-k');hold all
plot(h1,KUN(1,1),KUN(2,1),'ok')
plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
xlabel(h1,'\it\bfk')
ylabel(h1,'\it\bfU_N')

plot(h2,psis,'.-k');
set(h2,'Yscale','log','ylim',[.005 1])
ylabel('\it\bf\psi')
xlabel('\bfiterations')


%% Check starting point independence
figure(1302)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

alpha=1;
for I=1:100
    KUN=[0.003+0.005*sin(2*pi*I/100);0.005+0.015*cos(2*pi*I/100)]; %I.C.
    dKUN=[0.000004 0.00001]; %perturbations
    
    for III=1:10
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
        
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ)^-1*JJ'*psi2(III,1:11)';
    end
    
    plot(h1,KUN(1,:),KUN(2,:),'.-k');hold all
    plot(h1,KUN(1,1),KUN(2,1),'ok')
    plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,psis,'.-k');
    set(h2,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
    
end

%%
figure(1303)

set(gcf,'position',[60 60 680 660])
h1=axes('position',[.1 .05+1/3 .38 .27]);hold all
h2=axes('position',[.1 .05+2/3 .38 .27]);hold all
h4=axes('position',[.6 .05+1/3 .38 .27]);hold all
h3=axes('position',[.6 .05+2/3 .38 .27]);hold all
h5=axes('position',[.1 .05    .38 .27]);hold all

% contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

alpha=1;
KUN=[0.00; 0.0200; 10]; %I.C.
dKUN=[0.000004 0.00001 0.001]; %perturbations

for III=1:10
    CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,KUN(3,III),CC0); %model simulation
    psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
    psis(III)=norm(psi2(III,:));
    
    CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,KUN(3,III),CC0); %model simulation
    nCCtk=CCtk(TT+1)-CC;
    
    CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,KUN(3,III),CC0); %model simulation
    nCCtu=CCtu(TT+1)-CC;
    
    CCtv=CCtsim(KUN(1,III),KUN(2,III),t,UX,KUN(3,III)+dKUN(3),CC0); %model simulation
    nCCtv=CCtv(TT+1)-CC;
    
    JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1)  (nCCtu'-psi2(III,1:11)')/dKUN(2) (nCCtv'-psi2(III,1:11)')/dKUN(3)];
    
    KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ)^-1*JJ'*psi2(III,1:11)';
end
IIMv=[0.00206, 0.00221, 40.101];

plot(h1,KUN(1,:),KUN(2,:),'.-k');hold all
plot(h1,KUN(1,1),KUN(2,1),'ok')
plot(h1,IIMv(1),IIMv(2),'+b','markersize',10,'linewidth',2)
plot(h1,KUN(1,end),KUN(2,end),'xk','markersize',10,'linewidth',2)
xlabel(h1,'\it\bfk')
ylabel(h1,'\it\bfU_N')

plot(h2,KUN(1,:),KUN(3,:),'.-k');hold all
plot(h2,KUN(1,1),KUN(3,1),'ok')
plot(h2,IIMv(1),IIMv(3),'+b','markersize',10,'linewidth',2)
plot(h2,KUN(1,end),KUN(3,end),'xk','markersize',10,'linewidth',2)
xlabel(h2,'\it\bfk')
ylabel(h2,'\it\bfV')

plot(h3,KUN(3,:),KUN(2,:),'.-k');hold all
plot(h3,KUN(3,1),KUN(2,1),'ok')
plot(h3,IIMv(3),IIMv(2),'+b','markersize',10,'linewidth',2)
plot(h3,KUN(3,end),KUN(2,end),'xk','markersize',10,'linewidth',2)
xlabel(h3,'\it\bfV')
ylabel(h3,'\it\bfU_N')


plot3(h4,KUN(1,:),KUN(2,:),KUN(3,:),'.-k')
plot3(h4,KUN(1,1),KUN(2,1),KUN(3,1),'ok')
plot3(h4,IIMv(1),IIMv(2),IIMv(3),'+b','markersize',10,'linewidth',2)
plot3(h4,KUN(1,end),KUN(2,end),KUN(3,end),'xk','markersize',10,'linewidth',2)

view(h4,[35 35])
xlabel(h4,'\it\bfk')
ylabel(h4,'\it\bfU_N')
zlabel(h4,'\it\bfV')


plot(h5,psis,'.-k');
set(gca,'Yscale','log','ylim',[.005 1])
ylabel(h5,'\it\bf\psi')
xlabel(h5,'\bfiterations')


%% Check starting point independence
figure(1302)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

alpha=1;
for I=1:100
    KUN=[0.003+0.005*sin(2*pi*I/100);0.005+0.015*cos(2*pi*I/100)]; %I.C.
    dKUN=[0.000004 0.00001]; %perturbations
    
    for III=1:10
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
        
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ)^-1*JJ'*psi2(III,1:11)';
    end
    
    plot(h1,KUN(1,:),KUN(2,:),'.-k');hold all
    plot(h1,KUN(1,1),KUN(2,1),'ok')
    plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,psis,'.-k');
    set(h2,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
    
end

%%
figure(1305)
set(gcf,'position',[60 60 680 660])
h1=axes('position',[.1 .05+1/3 .38 .27]);hold all
h2=axes('position',[.1 .05+2/3 .38 .27]);hold all
h4=axes('position',[.6 .05+1/3 .38 .27]);hold all
h3=axes('position',[.6 .05+2/3 .38 .27]);hold all
h5=axes('position',[.1 .05    .38 .27]);hold all
% contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

alpha=1;

for I=1:100
            % these few lines set up a random constellation centred on
            % [.0025;.002;40] that sit on a crushed sphere shape
    KUN=.5-rand(3,1);
    KUN=KUN/norm(KUN).*[.0025;.002;20]+[.0025;.002;40];
    dKUN=[0.000004 0.00001 0.001]; %perturbations
    
    for III=1:10
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,KUN(3,III),CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,KUN(3,III),CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,KUN(3,III),CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        CCtv=CCtsim(KUN(1,III),KUN(2,III),t,UX,KUN(3,III)+dKUN(3),CC0); %model simulation
        nCCtv=CCtv(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1)  (nCCtu'-psi2(III,1:11)')/dKUN(2) (nCCtv'-psi2(III,1:11)')/dKUN(3)];
        
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ)^-1*JJ'*psi2(III,1:11)';
        KUN(1,III+1)=max(0,min(0.08,KUN(1,III+1)));
        KUN(2,III+1)=max(0,min(0.08,KUN(2,III+1)));
        KUN(3,III+1)=max(10,min(80,KUN(3,III+1)));
        
        
        
    end
    IIMv=[0.00206, 0.00221, 40.101];
    
    plot(h1,KUN(1,:),KUN(2,:),'.-k');hold all
    plot(h1,IIMv(1),IIMv(2),'+b','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xk','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,KUN(1,:),KUN(3,:),'.-k');hold all
    plot(h2,IIMv(1),IIMv(3),'+b','markersize',10,'linewidth',2)
    plot(h2,KUN(1,end),KUN(3,end),'xk','markersize',10,'linewidth',2)
    xlabel(h2,'\it\bfk')
    ylabel(h2,'\it\bfV')
    
    plot(h3,KUN(3,:),KUN(2,:),'.-k');hold all
    plot(h3,IIMv(3),IIMv(2),'+b','markersize',10,'linewidth',2)
    plot(h3,KUN(3,end),KUN(2,end),'xk','markersize',10,'linewidth',2)
    xlabel(h3,'\it\bfV')
    ylabel(h3,'\it\bfU_N')
    
    
    plot3(h4,KUN(1,:),KUN(2,:),KUN(3,:),'.-k')
    plot3(h4,IIMv(1),IIMv(2),IIMv(3),'+b','markersize',10,'linewidth',2)
    plot3(h4,KUN(1,end),KUN(2,end),KUN(3,end),'xk','markersize',10,'linewidth',2)
    
    view(h4,[35 35])
    xlabel(h4,'\it\bfk')
    ylabel(h4,'\it\bfU_N')
    zlabel(h4,'\it\bfV')
    
    plot(h5, psis,'.-k'); hold all
    set(gca,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
    
end

