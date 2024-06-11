clc
beep off
clear all
close all
format short g

set(0,'defaultfigureposition',[60 60 420,330])

t=0:600;                %mins
Un=0.01*ones(size(t));  %mmol/min
Ux=zeros(size(t));
Ux(62)=4;               %mmol
k=0.004;                %1/min
V=40;                   %L
C0=Un(1)/V/k;           %mmol/L

%analytical solution to the model
stime=zeros(1,3);
Ca=exp(-k*t).*(C0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));

%Picard (one possible numerical method) solution to the model
figure(501)
Cp=C0*ones(size(t));
plot(t,Ca,'k');hold all
for i=1:10  % a bunch of iterations in the Picard solution
                        
    plot(t,Cp,'r');     % a bunch of plotting lines
    title('\bfPicard')
    xlabel('\bftime [min]')
    ylabel('\bfCreatinine [mmolL^-^1]');
    axis([0 600 0 .2]);drawnow,pause(0.5)
    
    % the fucntional line:
    Cp=C0+cumtrapz(t,-k*Cp+(Ux+Un)/V);
    
end

% %Euler (another numerical method)
% Ce=zeros(size(t));
% Ce(1)=C0;
% dt=t(2)-t(1);
% hold off
% for tt=2:601
%     Ce(tt)=Ce(tt-1)*(1-k*dt)+(Un(tt-1)/V+Ux(tt-1)/V)*dt;
%     plot(t,Ca,'k',t,Cp,'r',t(1:tt),Ce(1:tt),'b');axis([0 600 0 .2]);
% xlabel('\bftime [min]')
% ylabel('\bfCreatinine [mmolL^-^1]')
% drawnow;%pause(0.0005)
% end
% legend('Analytical','Picard','Euler')

%% Into the real stuff:
figure(502)
set(gcf,'position',[60 60 420 330])
res=60;
                        
CC=Ca(1:res:601);       %data points
CC0=mean(CC(1:(60/res+1)));
TT=0:res:600;           %data times

CCt=interp1(TT,CC,t);   %plot data
plot(t/60,CCt,'g',TT/60,CC,'+k');hold all
axis([-.5 10.5 0 .20])

dCCt=diff(CCt);
dCCt(end+1)=dCCt(end);

figure(503)                                     
set(gcf,'position',[60 60 680 330]) % I want one this big
h1=axes('position',[.075 .1 .40 .85]);hold all  % set up some axis as handles so we can plot on differing figures as we wish
h2=axes('position',[.575 .1 .40 .85]);
figure(5031)
h3=axes;
% plot the current estimate:
plot(h3,[0 30],[0.004 0.004],[0 31],[0.01 0.01]);axis([0 31 0 .015]);hold all
A=[-CC' ones(600/res+1,1)/V ];
b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt(res+1:res:end)]'-[0 4/2/res 4/2/res zeros(1,600/res-2)]'/V;
%plot the cuurent mode simulation:
plot(h1,TT/60,CC,'+k');hold all;axis([-.5 10.5 0 .20]);drawnow

for I=1:30
    
    ss=A\b; % perform the inverse problem
            % resimulate the model at the new parameter values:   
    C=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(ss(2)+Ux)/V));
    
    psi(I)=norm(C(TT+1)-CC); %determine the residuals
    
    plot(h1,t/60,C,'r');hold all;xlabel(h1,'\bfTime [hours]');ylabel(h1,'\bfCreatinine [mmol/L]');
    plot(h2,psi,'o','markerfacecolor','b');xlabel(h2,'\bfiterations');ylabel(h2,'\bfresiduals \it\psi');drawnow
    plot(h3,I,ss(1),'.b',I,ss(2),'.g');xlabel(h3,'\bfiterations');ylabel(h3,'\bfvalues');hold all
            
    if I==1 % ask the observer for input before continuing
        input('press enter to continue')
    end
    
            % prepare for the next step
    dCCt=diff(C);
    dCCt(end+1)=dCCt(end);
        b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt(res+1:res:end)]'-[0 4/4 zeros(1,600/res-1)]'/V;
    
end
% put a legend on the relevant axis
legend(h3,'\itk_{true}','\itU_{N,true}','\itk_i','\itU_N_,_i')
%%
figure(504)

set(gcf,'position',[60 60 680 330])
h4=axes('position',[.075 .1 .40 .85]);hold all
h5=axes('position',[.575 .1 .40 .85]);
figure(5041)
h6=axes;
plot(h6,[0 30],[0.004 0.004],[0 31],[0.01 0.01]);axis([0 31 0 .015]);hold all

% to get the same outcomes as notes use this data:
CC=[0.0617362286206992,0.0578591828680488,0.138958525405719,0.126835540580120,0.121588113934300,0.105473590752212,0.0939101320790592,0.0905332846505024,0.0791227372065423,0.0721591272672233,0.0715330339740876;];
% else
% CC=Ca.*(1+0.02*randn(size(CC)));

CC0=mean(CC(1:(60/res+1)));
TT=0:res:600;

CCt=interp1(TT,CC,t);
plot(h4,t/60,CCt,'g',TT/60,CC,'+k');hold all
axis(h4,[-.5 10.5 0 .20])

dCCt=diff(CCt);
dCCt(end+1)=dCCt(end);

A=[-CC' ones(600/res+1,1)/V ];
b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt(res+1:res:end)]'-[0 4/2/60 4/2/60 zeros(1,600/res-2)]'/V;
            
for I=1:30  % use the same parameter ID methodology as above
    
    ss=A\b;
    C=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(ss(2)+Ux)/V));
    psi_n(I)=norm(C(TT+1)-CC); %determine the residuals
    
    plot(h4,t/60,C,'r');hold all;xlabel(h4,'\bfTime [hours]');ylabel(h4,'\bfCreatinine [mmol/L]');
    plot(h5,psi_n,'ob','markerfacecolor','b');xlabel(h5,'\bfiterations');ylabel(h5,'\bfresiduals \it\psi');drawnow
    plot(h6,I,ss(1),'.b',I,ss(2),'.g');xlabel(h6,'\bfiterations');ylabel(h6,'\bfvalues');hold all;hold all
  
    if I==1
        input('press enter to continue')
    end
    
    dCCt=diff(C);
    dCCt(end+1)=dCCt(end);
    b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt(res+1:res:end)]'-[0 4/4 zeros(1,600/res-1)]'/V;
    
end
legend(h6,'\itk_{true}','\itU_{N,true}','\itk_i','\itU_N_,_i')



%% Monte Carlo analysis
figure(505)
set(gcf,'position',[60 60 680 330])
h7=axes('position',[.085 .1 .40 .85]);hold all
h8=axes('position',[.585 .1 .40 .85]);

Ca=exp(-k*t).*(C0+cumtrapz(t,exp(k*t).*(Un+Ux)/V));
   
for jjj=1:50
            % add some random nosie to the perfect data
    CC=Ca(1:res:601).*(1+0.05*randn(1,600/res+1));
    CC0=mean(CC(1:(60/res+1)));
    TT=0:res:600;
    
    CCt=interp1(TT,CC,t);
    
    dCCt=diff(CCt);
    dCCt(end+1)=dCCt(end);
    
    A=[-CC' ones(600/res+1,1)/V ];
    b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt(res+1:res:end)]'-[0 4/2/60 4/2/60 zeros(1,600/res-2)]'/V;
                        
    for I=0:30          % do parameter ID
        
        ss=A\b;
        C=exp(-ss(1)*t).*(CC0+cumtrapz(t,exp(ss(1)*t).*(ss(2)+Ux)/V));
        psi(I+1)=norm(C(TT+1)-CC);
        dCCt=diff(C);
        dCCt(end+1)=dCCt(end);
        
        b=[dCCt(1) 0.5*dCCt(res:res:end)+0.5*dCCt(res+1:res:end)]'-[0 4/4 zeros(1,600/res-1)]'/V;
        
    end
    
    xlabel('\bf\itt \rm[hours]')
    ylabel('\bf\itC \rm[mmol/L]')
    axes(h7);plot(h7,psi,'.-b');hold all;axis([0 31 0 .04]);xlabel('\bfiterations');ylabel('\bfresiduals \it\psi')
    axes(h8);plot(h8,psi/psi(end),'.-b');hold all;axis([0 31 0.5 2]);xlabel('\bfiterations');ylabel('\bfnormalised residuals \it\psi/\psi')
   
end
