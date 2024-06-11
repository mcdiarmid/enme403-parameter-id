clc
beep off
clear all
close all
format short g
warning off
set(0,'units','centimeters','defaultfigureposition',[60 60 420 660])

t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62)=4;               %mmol
k=0.004;                %1/min
V=40;                   %L

TT=0:60:600;
%data with measurement error
CC=[0.0617362286206992,0.0578591828680488,0.138958525405719,0.126835540580120,0.121588113934300,0.105473590752212,0.0939101320790592,0.0905332846505024,0.0791227372065423,0.0721591272672233,0.0715330339740876;];
CC0=mean(CC(1:2));

%'perfect' forward simulation
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));

aa=linspace(0,0.008,31);                    %determine the grid resolution
bb=linspace(0,0.02,30);
psi=zeros(length(aa),length(bb));
xx=0;
for k_=aa
    xx=xx+1;yy=0;
    for UN_=bb
        yy=yy+1;
        CTT=fCCt(k_,UN_,t,Ux,CC0,V,TT+1);   %forward simulate at grid position
        psi(xx,yy)=norm(CTT-CC);            %store residual value
    end
end
figure(1701)
axes('position',[0.15 .58 .83 .4])
contour(aa,bb,psi',100);hold all            %plot objective surface
xlabel('\bf\itk')
ylabel('\bf\itU_N')

psis=psi*0;
x=[k,Un(1)]';       %I.C. 'parent' model values
for I=1:5           % perform Gauss Newton parameter ID
    plot(x(1),x(2),'+k')
    CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
    J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
    x=x-(J'*J)^-1*J'*(CT-CC');
    psis(I)=norm(CT-CC');
end
plot(x(1),x(2),'.k','markerfacecolor','w')

axes('position',[0.15 .08 .83 .4])          % position axis
plot(psis,'.-k')
xlabel('\bf\ititerations')
ylabel('\bf\it\psi')

xx=zeros(2,1000);
for II=1:length(xx)
    CC=C(TT+1).*(1+0.05*randn(size(TT)));   % get random noise on perfect simulation data
    x=[k,Un(1)]';
    for I=1:5                               % perform Gauss Newton parameter ID
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end

figure(1702)    % plot constellation of parameter values
set(gcf,'position',[60 60 420 330])
axis([0 .008 0 .02])
contour(aa,bb,psi',100);hold all
plot(xx(1,:),xx(2,:),'.k')
xlabel('\bf\itk')
ylabel('\bf\itU_N')
axis([0 .008 0 .02])


meank=mean(xx(1,:))     % find some summary statistics
CVk=std(xx(1,:))/meank
meanu=mean(xx(2,:))
CVu=std(xx(2,:))/meanu

[PCi,~,Contr]=pca(xx')  % principal component analysis
Contr=sqrt(Contr);      % sqrt it so it looks more like SD and less like variance!
100*Contr/sum(Contr)

figure(1703)
axes('position',[0.15 .58 .83 .4])
contour(aa,bb,psi',100);hold all
% plot mean SD outcomes
plot(meank*(1+CVk*randn(1000,1)),meanu*(1+CVu*randn(1000,1)),'.k')
axis([0 .008 0 .02])
xlabel('\bf\itk')
ylabel('\bf\itU_N')
rr1=randn(1000,1);
rr2=randn(1000,1);

axes('position',[0.15 .08 .83 .4])

contour(aa,bb,psi',100);hold all
% plot PCA outcomes
plot(meank+(Contr(1)*PCi(1)*rr1+Contr(2)*PCi(2,1)*rr2),meanu+(Contr(1)*PCi(1,2)*rr1+Contr(2)*PCi(2,2)*rr2),'.k')
axis([0 .008 0 .02])
xlabel('\bf\itk')
ylabel('\bf\itU_N')


%% start evaluating protocol alterations

% this section repeats a lot so I wont be re-annotating everything.
% better coding would have used functions for this, but ... um, shut up

% the original method
outcomes(1,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
figure(1706)
h1=plot(t,C,'k',TT,C(TT+1),'.k');hold all



figure(1704)
iC=cumtrapz(t,C);
mm(:,1:2)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];     % look at the A matrix formulation
conds(1)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);% calculate the condition number of the A matrix

set(gcf,'position',[60 60 680 660])
Ph1=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-k');hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-k')
xlabel('\bf\itk')
ylabel('\bf\itU_N')

figure(1707)
set(gcf,'position',[60 60 680 660])

subplot(3,3,1);
plot(t,t/(norm(TT)),'k',TT,TT/norm(TT),'.k');hold all
plot(t,iC/norm(iC(TT+1)),'k',TT,iC(TT+1)/norm(iC(TT+1)),'.k')

figure(1705)
set(gcf,'position',[60 60 680 660])

% the double data frequency protocol
t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62)=4;               %mmol
k=0.004;                %1/min
V=40;                   %L

TT=0:30:600;            % double the data frequency
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
iC=cumtrapz(t,C);
mm(1:21,3:4)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];  % look at the A matrix formulation
conds(2)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);% calculate the condition number of the A matrix

xx=zeros(2,1000);
for II=1:length(xx)     % same Monte Carlo method as above
    CC=C(TT+1).*(1+0.05*randn(size(TT)));
    x=[k,Un(1)]';
    for I=1:5
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end
% plot all the glorious outcomes!
figure(1706)
h2=plot(t,C,'b',TT,C(TT+1),'.b');hold all               % plot the data

figure(1707)
subplot(3,3,2);
plot(t,t/(norm(TT)),'b',TT,TT/norm(TT),'.b');hold all
plot(t,iC/norm(iC(TT+1)),'b',TT,iC(TT+1)/norm(iC(TT+1)),'.b')
figure(1705);
% plot the parameter constellation
subplot(3,2,1);plot(xx(1,:),xx(2,:),'.k');title('1');axis([0 .008 0 .02])
figure(1704);
[PCi,~,Contr]=pca(xx'); % plot the outcomes of the PCA analysis
Contr=sqrt(Contr);
outcomes(2,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
Ph2=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-b');hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-b')


% double the bolus
t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62)=8;               %mmol  %Double bolus
k=0.004;                %1/min
V=40;                   %L

TT=0:60:600;
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
iC=cumtrapz(t,C);
mm(1:11,5:6)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];
conds(3)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);

xx=zeros(2,1000);
for II=1:length(xx)
    CC=C(TT+1).*(1+0.05*randn(size(TT)));
    x=[k,Un(1)]';
    for I=1:5
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end

figure(1706)
h3=plot(t,C,'g',TT,C(TT+1),'.g');hold all 

figure(1707)
subplot(3,3,3);
plot(t,t/(norm(TT)),'g',TT,TT/norm(TT),'.g');hold all
plot(t,iC/norm(iC(TT+1)),'g',TT,iC(TT+1)/norm(iC(TT+1)),'.g')

figure(1705);
subplot(3,2,2);plot(xx(1,:),xx(2,:),'.k');title('2');axis([0 .008 0 .02])
figure(1704);
[PCi,~,Contr]=pca(xx');
Contr=sqrt(Contr);
outcomes(3,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
Ph3=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-g');hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-g')


% infusion rather than bolus
t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62:end)=0.02;        %mmol/min   % infusion
k=0.004;                %1/min
V=40;                   %L

TT=0:60:600;
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
iC=cumtrapz(t,C);
mm(1:11,7:8)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];
conds(4)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);

xx=zeros(2,1000);
for II=1:length(xx)
    CC=C(TT+1).*(1+0.05*randn(size(TT)));
    x=[k,Un(1)]';
    for I=1:5
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end
figure(1706)
h4=plot(t,C,'r',TT,C(TT+1),'.r');hold all

figure(1707)
subplot(3,3,4);
plot(t,t/(norm(TT)),'r',TT,TT/norm(TT),'.r');hold all
plot(t,iC/norm(iC(TT+1)),'r',TT,iC(TT+1)/norm(iC(TT+1)),'.r')
figure(1705);
subplot(3,2,3);plot(xx(1,:),xx(2,:),'.k');title('3');axis([0 .008 0 .02])
figure(1704);
[PCi,~,Contr]=pca(xx');
Contr=sqrt(Contr);
outcomes(4,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
Ph4=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-r');hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-r')


% Delayed bolus
t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(182)=4;              %mmol   %delay!
k=0.004;                %1/min
V=40;                   %L

TT=0:60:600;
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
iC=cumtrapz(t,C);
mm(1:11,9:10)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];
conds(5)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);

xx=zeros(2,1000);
for II=1:length(xx)
    CC=C(TT+1).*(1+0.05*randn(size(TT)));
    x=[k,Un(1)]';
    for I=1:5
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end
figure(1706)
h5=plot(t,C,'c',TT,C(TT+1),'.c');hold all

figure(1707)
subplot(3,3,5);
plot(t,t/(norm(TT)),'c',TT,TT/norm(TT),'.c');hold all
plot(t,iC/norm(iC(TT+1)),'c',TT,iC(TT+1)/norm(iC(TT+1)),'.c')
figure(1705);
subplot(3,2,4);plot(xx(1,:),xx(2,:),'.k');title('4');axis([0 .008 0 .02])
figure(1704);
[PCi,~,Contr]=pca(xx');
Contr=sqrt(Contr);
outcomes(5,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
Ph5=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-c');hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-c')


% half the measument error
t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62)=4;               %mmol
k=0.004;                %1/min
V=40;                   %L

TT=0:60:600;
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
iC=cumtrapz(t,C);
mm(1:11,11:12)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];
conds(6)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);

xx=zeros(2,1000);
for II=1:length(xx)
    CC=C(TT+1).*(1+0.025*randn(size(TT)));  % half mreasurement error
    x=[k,Un(1)]';
    for I=1:5
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end
figure(1706)
h6=plot(t,C,'-',TT,C(TT+1),'.','color',[1 .5 0]);hold all

figure(1707)
subplot(3,3,6);
plot(t,t/(norm(TT)),TT,TT/norm(TT),'.','color',[1 .5 0]);hold all
plot(t,iC/norm(iC(TT+1)),TT,iC(TT+1)/norm(iC(TT+1)),'.','color',[1 .5 0])

figure(1705);
subplot(3,2,5);plot(xx(1,:),xx(2,:),'.k');title('5');axis([0 .008 0 .02])
figure(1704);
[PCi,~,Contr]=pca(xx');
Contr=sqrt(Contr);
outcomes(6,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
Ph6=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-','color',[1 .5 0]);hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-','color',[1 .5 0])


% longer experiment
t=0:1200;               %mins   % longer
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62)=4;               %mmol
k=0.004;                %1/min
V=40;                   %L

TT=0:120:1200;          % longer
C=exp(-k*t).*(CC0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
iC=cumtrapz(t,C);
mm(1:11,13:14)=[iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')];
conds(7)=cond([iC(TT+1)'/norm(iC(TT+1)) TT'/norm(TT')]);

xx=zeros(2,1000);
for II=1:length(xx)
    CC=C(TT+1).*(1+0.05*randn(size(TT)));
    x=[k,Un(1)]';
    for I=1:5
        CT=fCCt(x(1),x(2),t,Ux,CC0,V,TT+1)';
        J=[fCCt(x(1)+0.000001,x(2),t,Ux,CC0,V,TT+1)'-CT fCCt(x(1),x(2)+0.000001,t,Ux,CC0,V,TT+1)'-CT]/.000001;
        x=x-(J'*J)^-1*J'*(CT-CC');
    end
    xx(:,II)=x;
end
figure(1706)
h7=plot(t,C,'m',TT,C(TT+1),'.m');hold all
% legend(h1(1:3:11),'.rig','.rig','2xfreq','2xfreq','2xUX','2xUX','infu','infu','delayUX','delayUX','.5noise',',5noise','20hrs','20hrs')
legend([h1(1) h2(1) h3(1) h4(1) h5(1) h6(1) h7(1)],'.rig','2xfreq','2xUX','infu','delayUX','.5noise','20hrs')
xlabel('\bfTime [mins]')
ylabel('\bfCreatinine [mmol.L^-^1]')
axis([0 1200 0 0.3])
set(gcf,'position',[60 60 680 660])

figure(1707)
subplot(3,3,7);
plot(t,t/(norm(TT)),'m',TT,TT/norm(TT),'.m');hold all
plot(t,iC/norm(iC(TT+1)),'m',TT,iC(TT+1)/norm(iC(TT+1)),'.m')

figure(1705);
subplot(3,2,6);plot(xx(1,:),xx(2,:),'.k');title('6');axis([0 .008 0 .02])
figure(1704);
[PCi,~,Contr]=pca(xx');
Contr=sqrt(Contr);
outcomes(7,:)=[mean(xx(1,:)) std(xx(1,:))/meank mean(xx(2,:)) std(xx(2,:))/meanu PCi(1,:) PCi(2,:) Contr'];
Ph7=plot(mean(xx(1,:))+[-1 0 1]*(PCi(1,1)*Contr(1)),mean(xx(2,:))+[-1 0 1]*(PCi(1,2)*Contr(1)),'.-m');hold all
plot(mean(xx(1,:))+[-1 0 1]*(PCi(2,1)*Contr(2)),mean(xx(2,:))+[-1 0 1]*(PCi(2,2)*Contr(2)),'.-m')

legend([Ph1(1) Ph2(1) Ph3(1) Ph4(1) Ph5(1) Ph6(1) Ph7(1)],'.rig','2xfreq','2xUX','infu','delayUX','.5noise','20hrs')


figure % this plots the norm of the matrix A columns. - Interesting from the point of view of the theory, but not really good for much else!
plot(mm(1:11,[1 5:2:end]));hold all
plot(mm(1:11,2),'k','linewidth',2)

