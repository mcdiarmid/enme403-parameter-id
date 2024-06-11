clc
beep off
clear all
close all
format short g

set(0,'defaultfigureposition',[60 60 420,150])
                        
tt=0:.05:100;           % define the temporal vectors
e1=zeros(size(tt));        e1(sort(floor(length(tt)*rand(1,10)+1)))=repmat([1 -1],1,5);e1=cumsum(e1);
e2=zeros(size(tt));        e2(sort(floor(length(tt)*rand(1,10)+1)))=repmat([1 -1],1,5);e2=cumsum(e2);
e3=zeros(size(tt));e3(1)=1;e3(sort(floor(length(tt)*rand(1,10)+1)))=repmat([-1 1],1,5);e3=cumsum(e3);
e4=e1.*e2;              % define the combination vectors
e5=e1.*e3;
e6=e2.*e3;
                                    
th=round(10*rand(6,1))/10;          % define the model values - round to nearest 0.1!

ee=[e1' e2' e3' e4' e5' e6'];
oc=ee*th;               % determine the outcomes
ocr=oc+randn(size(oc)); % add noise to the outcomes
th2=ee\ocr;             % do inverse problem
[th th2]
            
figure(401) % plot one of the functions            
area(tt,e1)
ylim([0 1.1])
xlabel('\it\bft')
ylabel('\it\bff_1')

figure(402) % plot the measured data
set(gcf,'position',[60 60 420,330])
plot(tt, ocr,'.k');hold all
xlabel('\it\bft')
ylabel('\it\bfX')

figure(403) % plot the identified outcomes with the parent model
set(gcf,'position',[60 60 420,330])
plot(tt, ocr,'.k');hold all
plot(tt, [ee*th ee*th2],'linewidth',2)
xlabel('\it\bft')
ylabel('\it\bfX')
legend('data','original','modelled')
