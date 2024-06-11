clear
clc
close all
beep off
warning off
format short g

set(0,'units','centimeters','defaultfigureposition',[60 60 680 330])

t = 0:.001:10;

L = 0.007;
R = 0.03;
Cinv = 1/1.4;  % 1/C

%generate some crazy voltage profile:
V = zeros(size(t));
V(1000:1000) = 30;
V(4000:5000) = .5;
V(5000:7000) = linspace(.5,-.5,2001);
V(7000:8000) = -.5;

figure(901);
subplot(1,2,1);
plot(t,V)
xlabel('\bfTime [s]')
ylabel('\bfVoltage across RLC [V]')
axis([0 10 -1 1])
set(gca,'ytick',-1:.5:1,'yticklabel',{'','-.5' '0' '.5','30'})

Vd = diff(V)/t(2);
t = t(1:end-1);

u = zeros(3,length(t));

%forward simulate the model: Simple Euler numerical integration. u(1,:) is
%the current in the circuit.
for I = 2:length(t)
    u(3,I) = Vd(I)/L -Cinv/L * u(1,I-1) - R/L * u(2,I-1); % d^2I/dt^2
    u(2,I) = u(2,I-1) + u(3,I)*t(2);  % dI/dt
    u(1,I) = u(1,I-1) + u(2,I)*t(2);  % I
end
subplot(1,2,2);
plot(t,u(1,:));hold all
xlabel('\bfTime [s]')
ylabel('\bfCurrent [A]')
figure(902)

% this is where we generate the matrix of interest:
% have a look at the matrix in the workspace viewer
A(:,1:3) = [[0 u(1,1:end-1)]' [0 0 u(1,1:end-2)]' [Vd(1:end-2) 0 0]'];
                        
spy(A(998:1005,:))      % spy plot where the entries of a sparse matrix are/
title('\bfPoints of matrix filled')
set(gca,'yticklabel',{' ' '998' '999','1000','1001','1002','1003','1004','1005',' '},'xticklabel',{'a1','a2','b'})

figure(903)
plot(t,u(1,:));hold all

b = [u(1,1:end-2) 0 0]';
[A(998:1005,:) b(998:1005)]

% limit identifiction to the 3 seconds of data (i.e. the first 3000 matrix inputs)
x = A(1:3000,:) \ b(1:3000)
uu = zeros(size(t));

% forward simulate the identified model
for I = 3:length(t)
    uu(I) = uu(I-1)*x(1) + uu(I-2)*x(2) + Vd(I)*x(3);
end

plot(t,uu);hold all
bb=A(1:3000,:)*x;       % the model for the period over which it was idntified
plot(t(1:length(bb)),bb)
legend('true','simulated','training')
xlabel('\bfTime [s]')
ylabel('\bfCurrent [A]')

%% ignore contributions from the model itself. this time the matrix A must be fully defined by the inputs.
AAA=zeros(3001,3000);
for I=1:3000    % this bit is a little counter intuitive, you should pay attention to the matrix in the workspace
    AAA(I:3001,I)=V(1:3002-I);
end
                                    
x=(AAA'*AAA)^-1*AAA'*b(1:3001);     % nope, it doesnt like it!!!
disp(['x=NaN? ' num2str(max(isnan(x)))]);
disp(['Fullrank? ' num2str(rank(AAA)==length(AAA))])

%get rid of zero rows:
bb=b(find(sum(abs(AAA),2)));
AA=AAA(find(sum(abs(AAA),2)),:);
AA=AA(:,find(sum(abs(AA))));

disp(['Fullrank? ' num2str(rank(AA)==length(AA))])
x=(AA'*AA)^-1*AA'*bb;               % great success
disp(['x=NaN? ' num2str(max(isnan(x)))]);

%% this reaction to the impulse can be used to simulate the response to all input signals
AAA=zeros(10001,length(x));
for I=1:length(x)
    AAA(I:10001,I)=V(1:10002-I);
end
uu=AAA*x;   % a sneaky approach to forward simulation
figure(904)
plot(t,u(1,:),[t 10],uu,t(1:3001),u(1,1:3001))
legend('true','simulated','training')
xlabel('\bfTime [s]')
ylabel('\bfCurrent [A]')

figure(905)
plot(t,u(1,:),[t 10],uu,t(1:3001),u(1,1:3001));axis([2 4 -.1 .1])
legend('true','simulated','training')
xlabel('\bfTime [s]')
ylabel('\bfCurrent [A]')

% try more data
AA=30*eye(3000);
bb=u(1,1000:3999)';
xx=(AA'*AA)^-1*AA'*bb;
AAA=zeros(10001,length(xx));
for I=1:length(xx)
    AAA(I:10001,I)=V(1:10002-I);
end
uu2=AAA*xx; % a sneaky approach to forward simulation

figure(906)
plot(t,u(1,:),[t 10],uu,[t 10],uu2,t(1:3999),u(1,1:3999))
legend('true','simulated 2s','simulated 3s','training')
xlabel('\bfTime [s]')
ylabel('\bfCurrent [A]')


