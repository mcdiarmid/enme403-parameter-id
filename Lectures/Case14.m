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
figure(1401)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

for lv=[0 .1 .3 1 3 10 30 100]
    
    alpha=1;
    
    KUN=[0.00; 0.0200]; %I.C.
    dKUN=[0.000004 0.00001]; %perturbations
    
    for III=1:30
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
            % Levenbergy goodness:
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ+lv*eye(2))^-1*JJ'*psi2(III,1:11)';
    end
    
    plot(h1,KUN(1,:),KUN(2,:),'.-','linewidth',2);hold all
    
    plot(h1,KUN(1,1),KUN(2,1),'ok')
    plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,psis,'.-');
    set(h2,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
end

legend(h2,{'\lambda=0' '0.1' '0.3' '1' '3' '10' '30' '100'})
figure(1402)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all

for lv=[0 .1 .3 1 3 10 30 100]
    
    alpha=1;
    
    KUN=[0.008; 0.000]; %I.C.
    dKUN=[0.000004 0.00001]; %perturbations
    
    for III=1:30
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
            % Levenbergy goodness
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ+lv*eye(2))^-1*JJ'*psi2(III,1:11)';
    end
    
    plot(h1,KUN(1,:),KUN(2,:),'.-','linewidth',2);hold all
    
    plot(h1,KUN(1,1),KUN(2,1),'ok')
    plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,psis,'.-');
    set(h2,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
end

legend(h2,{'\lambda=0' '0.1' '0.3' '1' '3' '10' '30' '100'})
%% Marquardt
figure(1403)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all
for lv=[0 .1 .3 1 3 10 30 100]
    
    alpha=1;
    KUN=[0.00; 0.0200]; %I.C.
    dKUN=[0.000004 0.00001]; %perturbations
    
    for III=1:30
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
            % Levenbergy-Marquardty Goodness
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ+lv*eye(2).*(JJ'*JJ))^-1*JJ'*psi2(III,1:11)';
    end
    
    plot(h1,KUN(1,:),KUN(2,:),'.-','linewidth',2);hold all
    plot(h1,KUN(1,1),KUN(2,1),'ok')
    plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,psis,'.-');
    set(h2,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
    
end
legend(h2,{'\lambda=0' '0.1' '0.3' '1' '3' '10' '30' '100'})

figure(1404)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.1 .12 .38 .85]);hold all
h2=axes('position',[.6 .12 .38 .85]);

contour(h1,-.002:0.0001:0.008,-.01:0.0002:0.02,log(psi'),50);hold all
for lv=[0 .1 .3 1 3 10 30 100]
                                    
    alpha=1;
    KUN=[0.008; 0.000];             %I.C.
    dKUN=[0.000004 0.00001];        %perturbations
    
    for III=1:30
        CCt=CCtsim(KUN(1,III),KUN(2,III),t,UX,V,CC0); %model simulation
        psi2(III,1:11)=CCt(TT+1)-CC;  % The model simulation (CCt) at the sample times (TT+1) minus the raw data
        psis(III)=norm(psi2(III,:));
        
        CCtk=CCtsim(KUN(1,III)+dKUN(1),KUN(2,III),t,UX,V,CC0); %model simulation
        nCCtk=CCtk(TT+1)-CC;
        
        CCtu=CCtsim(KUN(1,III),KUN(2,III)+dKUN(2),t,UX,V,CC0); %model simulation
        nCCtu=CCtu(TT+1)-CC;
        
        JJ=[(nCCtk'-psi2(III,1:11)')/dKUN(1) (nCCtu'-psi2(III,1:11)')/dKUN(2)];
            % Levenbergy-Marquardty Goodness
        KUN(:,III+1)=KUN(:,III)-alpha*(JJ'*JJ+lv*eye(2).*(JJ'*JJ))^-1*JJ'*psi2(III,1:11)';
    end
    
    plot(h1,KUN(1,:),KUN(2,:),'.-','linewidth',2);hold all
    plot(h1,KUN(1,1),KUN(2,1),'ok')
    plot(h1,0.00207, 0.002271,'+r','markersize',10,'linewidth',2)
    plot(h1,KUN(1,end),KUN(2,end),'xw','markersize',10,'linewidth',2)
    xlabel(h1,'\it\bfk')
    ylabel(h1,'\it\bfU_N')
    
    plot(h2,psis,'.-');
    set(h2,'Yscale','log','ylim',[.005 1])
    ylabel('\it\bf\psi')
    xlabel('\bfiterations')
    
end
legend(h2,{'\lambda=0' '0.1' '0.3' '1' '3' '10' '30' '100'})

