clc
beep off
clear all
close all
format short g

figure(601) % sizes and handles of axis/figures
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.085 .1 .40 .85]);hold all
h2=axes('position',[.585 .1 .40 .85]);

% parameter values
k=0.004;
V=40;
t=0:600;
UX=zeros(size(t));
UX(62)=4;
UN=0.01*ones(size(t));

res=60;
C=exp(-k*t).*(UN(1)/k/V+cumtrapz(t,exp(k*t).*(UN+UX)/V));
CC=C(1:res:end);
CC0=mean(CC(1:(60/res+1)));
TT=0:res:600;
CCt=interp1(TT,CC,t);

plot(h1,TT/60,CC,'+k');hold all;axis(h1,[0 10 0 .2])

% fill up the matrix with the coeeficients of the 3 parameter model
A=[-CC' ones(600/res+1,1) [0 4/2/res 4/2/res zeros(1,600/res-2)]'];
dCCt=diff(CCt);
dCCt(end+1)=dCCt(end);

for I=0:20  % iterate towards a solution
    
    b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt((res+1):res:end)]';
    ss=A\b;
    
    C=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(ss(2)+UX*ss(3))));
    psi(I+1)=norm(C(TT+1)-CC);
    plot(h1,t/60,C,'r');xlabel(h1,'\bf\itt \rm[hours]');ylabel(h1,'\bf\itC \rm[mmol/L]')
    plot(h2,psi,'o','markerfacecolor','b');xlabel('\bfiterations');ylabel('\bfresiduals \it\psi')
    
    if I==0
        input('press enter to continue')
    end
    dCCt=diff(C);
    dCCt(end+1)=dCCt(end);
    A=[-CC' ones(600/res+1,1) [0 1 zeros(1,600/res-1)]'];
    
end

figure(602)
set(gcf,'position',[60 60 680 330])
h3=axes('position',[.085 .1 .40 .85]);hold all
h4=axes('position',[.585 .1 .40 .85]);


k=0.004;
V=40;
t=0:600;
UX=zeros(size(t));
UX(62)=4;
UN=0.01*ones(size(t));

res=60;
C=exp(-k*t).*(UN(1)/k/V+cumtrapz(t,exp(k*t).*(UN+UX)/V));
CC=C(1:res:end);        % generate perfect data

for JJJ=1:100
            % add noise to data
    CC=C(1:res:end).*(1+0.05*randn(1,600/res+1));
    CC0=mean(CC(1:(60/res+1)));
    TT=0:res:600;
    CCt=interp1(TT,CC,t);
    
            % use the same parameter ID as before
    A=[-CC' ones(600/res+1,1) [0 4/2/res 4/2/res zeros(1,600/res-2)]'];
    dCCt=diff(CCt);
    dCCt(end+1)=dCCt(end);
    
    for I=0:30
        
        b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt((res+1):res:end)]';
        
        ss=A\b;
        ss(isnan(ss))=0;
        C=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(ss(2)+UX*ss(3))));
        
        psi(I+1)=norm(C(TT+1)-CC);
        
        dCCt=diff(C);
        dCCt(end+1)=dCCt(end);
        A=[-CC' ones(600/res+1,1) [0 1 zeros(1,600/res-1)]'];
        
    end
            %plot the convergence paths
    sss(JJJ,:)=[ss' psi(end)==min(psi)];
    axes(h3),plot(psi,'.-b','markerfacecolor','b');axis([0 31 0 .1]);xlabel('\bfiterations');ylabel('\bfresiduals \it\psi');hold all
    axes(h4),plot(psi/psi(end),'.-b','markerfacecolor','b');axis([0 31 0.5 2]);xlabel('\bfiterations');ylabel('\bfnormalised residuals \it\psi/\psi');hold all
end
sss(find(abs(sss)<0.00001))=0; %get rid of '-0'
%plot the rates of the parameter values
figure(603)
set(gcf,'position',[60 60 680 220])
h5=axes('position',[.09 .2 .22 .75]);hold all
h6=axes('position',[.09+1/3 .2 .22 .75]);hold all
h7=axes('position',[.09+2/3 .2 .22 .75]);hold all


plot(h5,sort(sss(:,1)),1:100);hold all
plot(h5,sort(sss(sss(:,4)==1,1)),1:length(sort(sss(sss(:,4)==1,1))));
xlim(h5,[-.01 .01]);
plot(h5,[k k],[0 100],'k');
ylabel(h5,'\bfrate (%)');
xlabel(h5,'\bf\it k')



plot(h6,sort(1./sss(:,3)),1:100);hold all
plot(h6,sort(1./sss(sss(:,4)==1,3)),1:length(sort(1./sss(sss(:,4)==1,3))));
xlim(h6,[-10000 10000]);
plot(h6,[V V],[0 100],'k');
ylabel(h6,'\bfrate (%)');
xlabel(h6,'\bf\it V')


plot(h7,sort(sss(:,2)./sss(:,3)),1:100);hold all
plot(h7,sort(sss(sss(:,4)==1,2)./sss(sss(:,4)==1,3)),1:length(sort(sss(sss(:,4)==1,2)./sss(sss(:,4)==1,3))));
xlim(h7,[-.1 .1])
plot(h7,[UN(1) UN(1)],[0 100],'k')
ylabel(h7,'\bfrate (%)');
xlabel(h7,'\bf\it U_N')

