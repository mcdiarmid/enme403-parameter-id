clc
beep off
clear all
close all
format short g
warning off

set(0,'defaultfigureposition',[60 60 420,330])
figure(701)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.085 .1 .40 .85]);hold all
h2=axes('position',[.585 .1 .40 .85]);
figure(702)
h3=axes;
plot(h3,[0 31],[0.004 0.004],'b',[0 31],[0.01 0.01],'g',[0 31],[1/40 1/40],'r');axis([0 11 0 .03]);hold all
                        
t=0:600;                % mins
Un=0.01*ones(size(t));  % mmol/min
Ux=zeros(size(t));
Ux(62)=4;               % mmol
k=0.004;                % 1/min
V=40;                   % L
C0=Un(1)/V/k;           % mmol/L

Ca=exp(-k*t).*(C0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));

TT=0:60:600;
CC=Ca(TT+1);
plot(h1,TT/60,CC,'+k');hold all


iUx=cumtrapz(t,Ux);     % integral of Ux
iUxi=iUx(TT+1);         % index of iUx
%TT is the integral of 1 w.r.t. time at the sample times!

Cg=zeros(size(t));      % bad guess of creatinine
b=CC'-mean(CC(1:2));


for III=1:10
    
    iCg=cumtrapz(t,Cg); % integrate the current Creatinine estimate
    iCgi=iCg(TT+1);     % get the values of the integral at the sample times
    A=[-iCgi' TT' iUxi'];           % fill up a matrix with the coefficient integrals at teh sample points
    
    ss=A\b;             
                        % forward simulate with new model values
    Cg=exp(-ss(1)*t).*(C0+cumtrapz(t,exp(ss(1)*t).*(Ux*ss(3)+ss(2))));
    
    psi(III)=norm(Cg(TT+1)-CC);     % store residuals
    plot(h1,t/60,Cg,'r');xlabel(h1,'\bf\itt \rm[hours]');ylabel(h1,'\bf\itC \rm[mmol/L]')
    plot(h2,psi,'ob','markerfacecolor','b');xlabel(h2,'\bfiterations');ylabel(h2,'\bfresiduals \it\psi')
    plot(h3,III,ss(1),'.b',III,ss(2)/ss(3),'.g',III,ss(3),'.r');xlabel(h3,'\bfiterations');ylabel(h3,'\bfparameter values')
    
    if III==1
        input('press enter to continue')
    end
end
legend(h3,'\itk','\itU_N','\itV^-^1')
th=[ss(1) ss(2)/ss(3) 1/ss(3)]      % print final parameter estimates

%% with noise
figure(703)
set(gcf,'position',[60 60 680 330])
h4=axes('position',[.085 .1 .40 .85]);hold all
h5=axes('position',[.585 .1 .40 .85]);
figure(704)
h6=axes;
plot(h6,[0 31],[0.004 0.004],'b',[0 31],[0.01 0.01],'g',[0 31],[1/40 1/40],'r');axis([0 11 0 .03]);hold all

CC=Ca(TT+1).*(1+0.05*randn(size(TT)));
% or to get parity with notes:
CC=[0.0617362286206992,0.0578591828680488,0.138958525405719,0.126835540580120,0.121588113934300,0.105473590752212,0.0939101320790592,0.0905332846505024,0.0791227372065423,0.0721591272672233,0.0715330339740876;];
CC0=mean(CC(1:2));

TT=0:60:600;
plot(h4,TT/60,CC,'+k');hold all

iUx=cumtrapz(t,Ux);     % integral of Ux
iUxi=iUx(TT+1);         % index of iUx
%TT is the integral of 1 w.r.t. at the sample times!

Cg=zeros(size(t));      % a bad guess of creatinine
b=CC'-CC0;

for III=1:10
    
    iCg=cumtrapz(t,Cg);
    iCgi=iCg(TT+1);
    A=[-iCgi' TT' iUxi'];
    
    ss=A\b;
    
    Cg=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(Ux*ss(3)+ss(2))));
    
    psi(III)=norm(Cg(TT+1)-CC);
    plot(h4,t/60,Cg,'r');xlabel(h4,'\bf\itt \rm[hours]');ylabel(h4,'\bf\itC \rm[mmol/L]')
    plot(h5,psi,'ob','markerfacecolor','b');xlabel(h5,'\bfiterations');ylabel(h5,'\bfresiduals \it\psi')
    plot(h6,III,ss(1),'.b',III,ss(2)/ss(3),'.g',III,ss(3),'.r');xlabel(h6,'\bfiterations');ylabel(h6,'\bfparameter values')
    
    if III==1
        input('press enter to continue')
    end
end
legend(h6,'\itk','\itU_N','\itV^-^1')
th=[ss(1) ss(2)/ss(3) 1/ss(3)]

%% Monte Carlo

figure(705)
set(gcf,'position',[60 60 680 330])
h7=axes('position',[.085 .1 .40 .85]);hold all
h8=axes('position',[.585 .1 .40 .85]);
h81=axes('position',[.80 .62 .18 .35]);hold all


k=0.004;
V=40;
t=0:600;
TT=0:60:600;
UX=zeros(size(t));
UX(62)=4;
UN=0.01*ones(size(t));

C=exp(-k*t).*(UN(1)/k/V+cumtrapz(t,exp(k*t).*(UN+UX)/V));

for JJJ=1:100
                                                
    CC=C(TT+1).*(1+0.05*randn(size(TT)));       % add noise to the data
    CC0=mean(CC(1:2));
    Cg=zeros(size(t));  % start form the same vector everytime
    b=CC'-CC0;
                        
    for III=1:10
        
        iCg=cumtrapz(t,Cg);
        iCgi=iCg(TT+1);
        A=[-iCgi' TT' iUxi'];
        
        ss=A\b;
        Cg=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(Ux*ss(3)+ss(2))));
        psi(III)=norm(Cg(TT+1)-CC);
    end
            % plot convergence paths
    sss(JJJ,:)=[ss' (psi(end)*.995)<min(psi)];
    axes(h7),plot(psi,'.-b','markerfacecolor','b');axis([0 11 0 .1]);xlabel('\bfiterations');ylabel('\bfresiduals \it\psi');hold all
    axes(h8),plot(psi/psi(end),'.-b','markerfacecolor','b');axis([0 11 0.5 2]);xlabel('\bfiterations');ylabel('\bfnormalised residuals \it\psi/\psi');hold all
    axes(h81),plot(psi/psi(end),'.-b','markerfacecolor','b');axis([0 11 0.995 1.005]);

end
figure(706)
set(gcf,'position',[60 60 680 220])
h9=axes('position',[.09 .2 .22 .75]);hold all
h10=axes('position',[.09+1/3 .2 .22 .75]);hold all
h11=axes('position',[.09+2/3 .2 .22 .75]);hold all


plot(h9,sort(sss(:,1)),1:100);hold all
plot(h9,sort(sss(sss(:,4)==1,1)),1:length(sort(sss(sss(:,4)==1,1))));
xlim(h9,[0 .01]);ylim(h9,[0 100]);
plot(h9,[k k],[0 101],'k');
ylabel(h9,'\bfrate (%)');
xlabel(h9,'\bf\it k')


plot(h10,sort(1./sss(:,3)),1:100);hold all
plot(h10,sort(1./sss(sss(:,4)==1,3)),1:length(sort(1./sss(sss(:,4)==1,3))));
xlim(h10,[30 50]);ylim(h10,[0 100]);
plot(h10,[V V],[0 101],'k');
ylabel(h10,'\bfrate (%)');
xlabel(h10,'\bf\it V')


plot(h11,sort(sss(:,2)./sss(:,3)),1:100);hold all
plot(h11,sort(sss(sss(:,4)==1,2)./sss(sss(:,4)==1,3)),1:length(sort(sss(sss(:,4)==1,2)./sss(sss(:,4)==1,3))));
xlim(h11,[0 .02]);ylim(h11,[0 100]);
plot(h11,[UN(1) UN(1)],[0 101],'k')
ylabel(h11,'\bfrate (%)');
xlabel(h11,'\bf\it U_N')

%% shit is about to get real!
set(0,'units','centimeters','defaultfigureposition',[60 60 650 650])
figure(707)
plot3(sss(:,1),sss(:,2),sss(:,3),'+k');hold all
nnn=100;
[x,y,z]=sphere(nnn); %define 100 points of a sphere
x=x/2;
y=y/2;
z=z/2;
x=(x+1)*0.004; %scale the sphere to reasonable bounds on the parameter space
y=(y+1)*0.00025;
z=(z+1)*0.025;

% find the residual error of the model when the points on the sphere are
% used in the model
for ggg=1:nnn+1 
    for hhh=1:nnn+1
        CCt=exp(-x(ggg,hhh)*t).*(CC0+cumtrapz(t,exp(x(ggg,hhh)*t).*(y(ggg,hhh)+UX*z(ggg,hhh))));
         err(ggg,hhh)=norm(CCt(TT+1)-CC);
    end
end

subplot(2,2,2);                                    
plot3(sss(:,1),sss(:,2),sss(:,3),'+k');hold all             % plot the parameter constellations
surf(x,y,z,log(err),'linestyle','none');alpha(0.5)          % plot the error sphere
xlabel('\bf\itk');ylabel('\bf\itU_N/V');zlabel('\bf\it1/V');
axis([0 0.008 0 0.0006 0 0.04])
view([0 0]) % choose some angle to view it at
            
subplot(2,2,3);
plot3(sss(:,1),sss(:,2),sss(:,3),'+k');hold all
surf(x,y,z,log(err),'linestyle','none');alpha(0.5)
xlabel('\bf\itk');ylabel('\bf\itU_N/V');zlabel('\bf\it1/V');
view([0 90]);axis([0 0.008 0 0.0006 0 0.04])


subplot(2,2,4);
plot3(sss(:,1),sss(:,2),sss(:,3),'+k');hold all
surf(x,y,z,log(err),'linestyle','none');alpha(0.5)
xlabel('\bf\itk');ylabel('\bf\itU_N/V');zlabel('\bf\it1/V');
view([90 0])
axis([0 0.008 0 0.0006 0 0.04])

subplot(2,2,1);
plot3(sss(:,1),sss(:,2),sss(:,3),'+k');hold all
surf(x,y,z,log(err),'linestyle','none');alpha(0.5)
xlabel('\bf\itk');ylabel('\bf\itU_N/V');zlabel('\bf\it1/V');
view([-22 30])
axis([0 0.008 0 0.0006 0 0.04])

%principal component analysis
[bb,cc,dd]=pca(sss(:,1:3));
var_prop=dd/sum(dd);
[bb' var_prop*100]
