clc
beep off
clear all
close all
format short g

set(0,'units','centimeters','defaultfigureposition',[60 60 680 330])
load PatientDataSet.mat % load in the patient data set (flow, pressure, time)

paw=paw_SA'; %rename everytime to my preferences!
Vdot=flow_SA';
time=t_SA;
V=cumtrapz(time,Vdot); %integrate flow for volume

figure(801)
axes('position',[.08 .12 .40 .85])
plot(time,paw,'k');hold all %plot the raw pressure data
axis([0 6 -1 20])
ylabel('\bf\itp_a_w \rm[cmH_20]')
xlabel('\bfTime \rm[s]')
axes('position',[.58 .12 .40 .85])
plot(time,V,'k');hold all %plot the raw volume data
axis([0 6 0 500])
ylabel('\bf\itV \rm[mL]')
xlabel('\bfTime \rm[s]')

%% FOM
A=[V' Vdot'];
b=paw';
ss=(A'*A)^-1*A'*b;
pawFOM=A*ss;
figure(802)
% plot(time,[paw;pawFOM'])
disp(ss)

axes('position',[.08 .12 .40 .85])
plot(time,paw,'k');hold all %plot the raw pressure data
plot(time,pawFOM,'r'); %plot the raw pressure data
axis([0 6 -1 20])
ylabel('\bf\itp_a_w \rm[cmH_20]')
xlabel('\bfTime \rm[s]')
axis([0 6 0 20])
axes('position',[.58 .12 .40 .85])
plot(time,zeros(size(time)),'k');hold all %plot the raw volume data
plot(time,pawFOM-paw','r');hold all %plot the raw volume data
axis([0 6 -5 5])
ylabel('\bf residuals\rm[cmH_20]')
xlabel('\bfTime \rm[s]')

%%
figure(803)
xx=0;yy=0;  % construct a residual surface
for E=linspace(0,0.05,50)
    xx=xx+1;yy=0;
    for R=linspace(0,0.02,51);      % always set differing vector lengths!
        yy=yy+1;
        err(xx,yy)=norm(A*[E;R]-b);
    end
end

set(gcf,'position',[60 60 420 330])
contour(linspace(0,0.05,50),linspace(0,0.02,51),err');hold all
plot(ss(1),ss(2),'+k')  % plot the identified parameter set against the model ojective surface
xlabel('\bfElastance [cmH_20/L]')
ylabel('\bfResistance [cmH_20.s/L]')
                        

%% Viscoelastic
figure (804)
RHS=paw;
subplot(1,2,1);plot(time,paw,'k');hold all
for I=1:20
    a=Vdot;
    b=cumtrapz(time,paw);
    c=V;
    d=cumtrapz(time,V);
                                                        
    solns=lsqlin([a(:) b(:) c(:) d(:)],RHS(:));     % this is another way to do A\b etc. it offers options such as limits
    
    R1_=solns(1);
    r2c2=-1/(solns(2));
    C1_=1/(solns(4)*r2c2);
    C2_=1/(solns(3)-1/C1_-R1_/r2c2);
    R2_=r2c2/C2_;
    % the following line can be used to simulate the model in an explicit way.
    % It may be that a implicit calculation would be better, more realistic. In
    % such cases the commented out text below can be ressurected. note that
    % in this case there is virtually no difference between the implicit and explicit
    % simulations due to the limited changes in paw. 
    paw=solns(1)*a+solns(2)*b+solns(3)*c+solns(4)*d;
    
    
    %     p1dot_=Vdot/C1_;
    %     p1=cumtrapz(time,p1dot_);
    %     p2=exp(cumtrapz(time,-ones(1,length(time))/(R2_*C2_))).*cumtrapz(time,exp(cumtrapz(time,ones(1,length(time))/(R2_*C2_))).*Vdot/C2_);
    %
    %     for J=1:50
    %         p2=-1/R2_/C2_*cumtrapz(time,p2)+V/C2_;
    %     end
    %
    %     paw=p1+p2+R1_*Vdot;
    
    AAA(I,:)=([R1_ R2_ C1_ C2_ norm(paw-RHS)]);     % store the values
    
end
% axes('position',[a b c d]) is much better than subplot, nonetheless...
subplot(1,2,1);plot(time,paw,':b');hold all;ylim([0 20])
p1dot_=Vdot/C1_;        % evaluate the model
p1=cumtrapz(time,p1dot_);
p2=exp(cumtrapz(time,-ones(1,length(time))/(R2_*C2_))).*cumtrapz(time,exp(cumtrapz(time,ones(1,length(time))/(R2_*C2_))).*Vdot/C2_);

for J=1:50
    p2=-1/R2_/C2_*cumtrapz(time,p2)+V/C2_;
end

paw=p1+p2+R1_*Vdot;

AAA(I,:)=([R1_ R2_ C1_ C2_ norm(paw-RHS)]);
disp(AAA)
% plot all of the outcomes
subplot(1,2,1);plot(time,paw,'b');hold all
ylabel('\bf\itp_a_w \rm[cmH_20]')
xlabel('\bfTime \rm[s]')
set(gcf,'position',[60 60 680 330])

subplot(5,2,2);plot(AAA(:,1),'.-b');ylabel('\bf\itR_1')
subplot(5,2,4);plot(AAA(:,2),'.-b');ylabel('\bf\itR_2')
subplot(5,2,6);plot(AAA(:,3),'.-b');ylabel('\bf\itC_1')
subplot(5,2,8);plot(AAA(:,4),'.-b');ylabel('\bf\itC_2')
subplot(5,2,10);plot(AAA(:,5),'.-b');ylabel('\bf\it\psi');xlabel('\bfIterations');ylim([8.5 8.8])
            
figure(805) % plot the residuals
h7=axes('position',[.085 .12 .40 .85]);hold all

plot(time,pawFOM','r');hold all
plot(time,paw,'b');
plot(time,RHS,'k');
ylim([0 20])
ylabel('\bf\itp_a_w \rm[cmH_20]')
xlabel('\bfTime \rm[s]')

h8=axes('position',[.585 .12 .40 .85]);hold all

plot(time,pawFOM'-RHS,'r');hold all
plot(time,paw-RHS,'b');
plot(time,RHS-RHS,'k');
ylabel('\bf\itp_a_w \rm\bfresiduals \rm[cmH_20]')
xlabel('\bfTime \rm[s]')
plot([0 6],[0 0],'k');hold all
axis([0 6 -4 4])
legend('FOM','Viscoelastic')
