clc
beep off
clear all
close all
format short g

set(0,'units','centimeters','defaultfigureposition',[60 60 630 350])

LT=2000; %determine the simulation time
t=0:LT;
%% this is the hidden model, it will be largely ignored by the identification system
k=0.1+0.9*rand(8); % define a matrix of random values
k=k.*[  0 0 0 0 0 0 0 0; % eliminate teh elements of the matrix that do not represent a transport
        0 0 0 1 0 0 0 0;
        0 1 0 0 0 0 0 0;
        1 1 0 0 0 0 1 0;
        0 0 1 1 0 0 0 0;
        0 0 0 0 0 0 0 0;
        0 0 0 1 1 1 0 1;
        0 0 0 0 1 0 1 0];
k=k-diag(sum(k,1)+rand(1,8).*[0 1 0 0 1 0 1 1]); %add the outflow of each compartment, and allow 1st order clearance from comp 2, 5, 7 and 8

k=round(1000*k)/10000; %get nice roundish numbers

%% this is the inputs - each more darstardly than the last!
                                                            
I=zeros(8,LT+1);                                           
I(1,round(LT*rand(10,1)))=100;                              % determine a bunch of impulse measurements
I(2,:)=interp1(linspace(0,LT,11),rand(11,1),0:LT);          % determine random linear behaviour
I(3,:)=rand(1)*(1+sin(rand(1)*linspace(0,20,LT+1)+pi*rand(1)));         % some sine curve
ind=round(LT*sort(rand(10,1)));
I(6,ind)=(repmat([1 -1],1,5));
I(6,:)=cumsum(I(6,:));                                      % some random box function

figure(1001);
% set(gcf,'position',[60 60 420 350])
subplot(4,1,1);plot(t,I(1,:));ylim([0 105]);ylabel('\bfInput #1')
subplot(4,1,2);plot(t,I(2,:));ylabel('\bfInput #2')
subplot(4,1,3);plot(t,I(3,:));ylabel('\bfInput #3')
subplot(4,1,4);plot(t,I(6,:));ylim([0 1.1]);ylabel('\bfInput #4')
xlabel('\bfTime [s]')
% you can use these 'training-wheel' inputs if you like
% I=zeros(8,LT+1);
% I(1,10)=1;
% I(2,210)=1;
% I(3,610)=1;
% I(6,410)=1;

%% forward simulate the model using good old (bad) Euler method!
u=zeros(8,LT+1);
for ttt=2:LT+1
 u(:,ttt)=u(:,ttt-1)+k*u(:,ttt-1)+I(:,ttt);
end

figure(1002);
plot(0:LT,u(7:8,:)')
xlabel('\bfTime [s]')
ylabel('\bfOutputs');
legend('output #1','output #2')

%% let the parameter ID begin!!!
A=zeros(LT+1,1000);
                        
for jj=1:250            % determine a matrix the recors 250 seconds of input lag for each input over the full length of time
    A(jj:end,jj)=I(1,1:end-jj+1)';
    A(jj:end,250+jj)=I(2,1:end-jj+1)';
    A(jj:end,500+jj)=I(3,1:end-jj+1)';
    A(jj:end,750+jj)=I(6,1:end-jj+1)';
end
figure(1003)
spy(A)% have a look at the rank of the matrix

b7=u(7,:)';
b8=u(8,:)';

ID_p=1200; % detemine how much input data you want to use 1 second to LT seconds, <100 will fail everytime!!

ss7=A(1:ID_p,:)\b7(1:ID_p); % we don't want to use (A'A)^-1*A'*b because it goes unstable at such sensitivities
ss8=A(1:ID_p,:)\b8(1:ID_p);
% Note that we used the same input matrix for both analyses of this model!


%% simulate for all time using our model that was identified
uu=zeros(2001,2);

for ttt=2:LT
    uu(ttt,1)=A(ttt,:)*ss7;
    uu(ttt,2)=A(ttt,:)*ss8;
end
figure(1004);hold all
plot(t,u(7:8,:))
plot(0:LT,uu',':')
xlabel('\bfTime [s]')
ylabel('\bfOutputs');
legend('1true','2true','1recov','2recov')

figure(1005)
plot(t,u)
xlabel('\bfTime [s]')
ylabel('\bfOutputs');
legend('1','2','3','4','5','6','7 O#1','8 O#2')
