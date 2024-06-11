%% This is the first in a series of codes that are intended to show how the
% figures in the parameter identification notes were derived. You don't
% need to study or replicate these particular codes, but understanding them
% will contribute to your coding endevours. I will be the first to admit
% that the codes are sparsely annotated, and in places arent annotated at
% all. I have tried to annotate the key points that might form a barrier in
% understanding, but not much else. If you dont recognise a function try
% F1! If youre still having problems, unsupress the inputs and outputs.
% I also tend to write long scripts rather than a series of functions. This
% is not good coding practice. I can get away with it because the codes are
% relatively simple and would only be partially improved with functions. It
% also helps people to follow what is going on.
% Ultimately, if you cant understand what is going on, please print and
% bring the codes to lectures or to my office for explanations!

clc
beep off            % I hate the matlab beep that stings with every error
clear all
close all
format short g      % make numbers pretty

% get a better default figure
set(0,'units','centimeters','defaultfigureposition',[60 60 420,330])

dt=0.1;
tt=0:dt:20;         % some vector of time
xx=2.2*sqrt(tt)+2;  % the true model

figure(101)
X=[3.56 4.26 5.10 5.76 6.91 7.17 9.00]; % get some 'random' data
T=[0.5  1    2    3.2  4.8  6    10];   % sample times

% unsurpress for truely random data
% X=round(100*(xx(T*10+1).*(1+0.02*randn(1,7))))/100

plot(T,X,'+k');hold all
xlabel('\it\bft')
ylabel('\it\bfX')
axis([0 20 0 11])

ss=polyfit(T,X,6)'; % use matlabs function to find the coefficients of 6th order polyn
disp(['6th O polyn vals: ' num2str(ss')])
plot(tt,polyval(ss,tt))
disp(['6th O polyn @ t=8: ' num2str(polyval(ss,8))])

plot(8,polyval(ss,8),'or');     % evaluate the model at t=8s
xlabel('\it\bft')
ylabel('\it\bfX')
axis([0 12 0 12])

%%
figure(102)
ss=polyfit(T,X,1)';             % evaluate the 1st order polyn
disp(' ')
disp(['1st O polyn vals: ' num2str(ss')])
plot(T,X,'+k',tt,polyval(ss,tt),8,polyval(ss,8),'or');hold all
disp(['1st O polyn @ t=8: ' num2str(polyval(ss,8))])
xlabel('\it\bft')
ylabel('\it\bfX')
axis([0 12 0 12])

%% D.I.Y.
figure(103)
set(gcf,'position',[60 60 420,660])
A=[T' [1 1 1 1 1 1 1]'];
b=X';
disp('A=')
disp(A)
disp('b=')
disp(b)
ss2=(A'*A)^-1*A'*b;

axes('position',[.11 .56 .85 .42])
plot(T,X,'+k',tt,polyval(ss2,tt),'-r');hold all     % plot the model
for i=1:length(T)               % plot the residuals
    plot([T(i) T(i)],[X(i); polyval(ss2,T(i))],'-r');hold all
end
xlabel('\it\bft')
ylabel('\it\bfX')
axis([0 12 0 12])

axes('position',[.11 .06 .85 .42])
for i=1:length(T)               % plot the residuals again
    plot([T(i) T(i)],[0; X(i)-polyval(ss2,T(i))],'-xr');hold all
end
plot([0 12],[0 0],'color',[.8 .8 .8])               % plot the zero line
xlabel('\it\bft')
ylabel('\it\bf\psi')
axis([0 12 -1 1])

disp(['residual error= ' num2str(norm(X-polyval(ss2,T)))])
%%
figure(104)
set(gcf,'position',[60 60 420,660])

A=[sqrt(T') [1 1 1 1 1 1 1]'];
b=X';
ss3=(A'*A)^-1*A'*b;

disp(['sqrt fcn coeffs: ' num2str(ss3')])
disp(['residual error= ' num2str(norm(A*ss3-b))])

axes('position',[.11 .56 .85 .42])
plot(T,X,'+k',tt,ss2(1)*tt+ss(2),'r');hold all      % plot the data and 1st order model
plot(tt,ss3(1)*sqrt(tt)+ss3(2),'color',[0 .5 0])    % plot the sqrt model
set(0,'units','centimeters','defaultfigureposition',[60 60 750,330])
axis([0 12 0 12])
xlabel('\it\bft')
ylabel('\it\bfX')

axes('position',[.11 .06 .85 .42])
for i=1:length(T)               % plot and compare model residuals
    plot([T(i) T(i)],[0; X(i)-polyval(ss2,T(i))],'-xr');hold all
    plot([T(i) T(i)]+.1,[0; X(i)-ss3(1)*T(i)^.5-ss3(2)],'-x','color',[0 .5 0]);hold all
end
plot([0 12],[0 0],'color',[.8 .8 .8])
xlabel('\it\bft')
ylabel('\it\bf\psi')
axis([0 12 -1 1])

%%
figure (105)

for ss31=0:0.1:5                % define posiitons on grid
    for ss32=0:.1:5
        XXX=ss31*sqrt(0:0.1:10)+ss32;               % determine the model at points on the grid
        % determine the residual error at the points on the grid
        psi(round(ss31*10+1),round(ss32*10+1))=norm(XXX(T*10+1)-X);
    end
end
set(gcf,'position',[60 60 840,330])
axes('position',[.06 .15 .42 .82])
[cc,hh]=contourm(0:.1:5,0:.1:5,psi',1:2:21);hold all% plot the contours
clegendm(cc,hh,1)               % adds a legend to contour plots
plot(ss3(1),ss3(2),'k+','markersize',10,'linewidth',2)
xlabel('\it\bf\theta_1')
ylabel('\it\bf\theta_2')

axes('position',[.56 .15 .42 .82])
surf(0:.1:5,0:.1:5,psi','linestyle','none')         % plot nice 3D objective valley
xlabel('\it\bf\theta_1')
ylabel('\it\bf\theta_2')
zlabel('\it\bf\psi')
view([-6 44])
set(0,'units','centimeters','defaultfigureposition',[60 60 420,330])


%% try 2nd and 3rd order model
A=[T'.^3 T'.^2 T' ones(7,1)];   % 3rd order
b=X';
ss4=(A'*A)^-1*A'*b;
norm(A*ss4-b);

disp(['3rd order fcn coeffs: ' num2str(ss4')])
disp(['residual error= ' num2str(norm(A*ss4-b))])
A=[ T'.^2 T' ones(7,1)];        % 2nd order
b=X';
ss4=(A'*A)^-1*A'*b;
norm(A*ss4-b);

disp(['2nd order fcn coeffs: ' num2str(ss4')])
disp(['residual error= ' num2str(norm(A*ss4-b))])


%% MC analysis
figure(106)
x=2.2*sqrt(tt)+2; % the true model
X=xx(round(T/dt+1)); 

% detemine the coefficient matrices
A1=[T' [1 1 1 1 1 1 1]'];
A2=[T.^2' T' [1 1 1 1 1 1 1]'];
A3=[T.^3' T.^2' T' [1 1 1 1 1 1 1]'];
A4=[sqrt(T') [1 1 1 1 1 1 1]'];

% predefine some matrices to hold outcomes
RR1=zeros(10000, 8);
RR2=zeros(10000, 8);
RR3=zeros(10000, 8);
RR4=zeros(10000, 8);

for III=1:10000
    XX=X.*(1+0.02*randn(size(X)));  % add noise to the clincial data
    b=XX';
    ss1=A1\b;
    ss2=A2\b;
    ss3=A3\b;
    ss4=A4\b;
    RR1(III,:)=[A1*ss1-XX'; norm(A1*ss1-XX')];  % store the model residuals
    RR2(III,:)=[A2*ss2-XX'; norm(A2*ss2-XX')];
    RR3(III,:)=[A3*ss3-XX'; norm(A3*ss3-XX')];
    RR4(III,:)=[A4*ss4-XX'; norm(A4*ss4-XX')];
    
end

% perform summary statistics on the residuals - proper analysis would use
% standard error not standard deviation. Standard error is equal to
% standard deviation divided by the squareroot of the number of samples.
% This gives an idea of the certainty of the mean values, rather than the
% degree of variance.
RRR1=-[mean(RR1)-std(RR1) ;mean(RR1);mean(RR1)+std(RR1) ];
RRR2=-[mean(RR2)-std(RR2) ;mean(RR2);mean(RR2)+std(RR2) ];
RRR3=-[mean(RR3)-std(RR3) ;mean(RR3);mean(RR3)+std(RR3) ];
RRR4=-[mean(RR4)-std(RR4) ;mean(RR4);mean(RR4)+std(RR4) ];

errorbar(T,RRR1(2,1:7),RRR1(2,1:7)-RRR1(1,1:7),RRR1(2,1:7)-RRR1(3,1:7),'b');hold all
errorbar(T+.1,RRR2(2,1:7),RRR2(2,1:7)-RRR2(1,1:7),RRR2(2,1:7)-RRR2(3,1:7),'r');hold all
errorbar(T+.2,RRR3(2,1:7),RRR3(2,1:7)-RRR3(1,1:7),RRR3(2,1:7)-RRR3(3,1:7),'color',[0 .5 0]);hold all
errorbar(T+.3,RRR4(2,1:7),RRR4(2,1:7)-RRR4(1,1:7),RRR4(2,1:7)-RRR4(3,1:7),'k');hold all

MAphi= [mean(mean(abs(RR1(:,1:7)))) mean(mean(abs(RR2(:,1:7)))) mean(mean(abs(RR3(:,1:7))) ) mean(mean(abs(RR4(:,1:7))))]
MAMphi=[mean(abs(mean(RR1(:,1:7)))) mean(abs(mean(RR2(:,1:7)))) mean(abs(mean(RR3(:,1:7)))) mean(abs(mean(RR4(:,1:7))))]

legend('1st order','2nd order','3rd order','sqrt model')
axis([0 12 -1 1])
xlabel('\it\bft')
ylabel('\it\bf\psi')








