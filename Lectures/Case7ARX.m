clc
beep off
clear all
close all
format short g
warning off

set(0,'defaultfigureposition',[60 60 420,330])
figure(711)
% set(gcf,'position',[60 60 680 330])

CC=[0.0617362286206992,0.0578591828680488,0.138958525405719,0.126835540580120,0.121588113934300,0.105473590752212,0.0939101320790592,0.0905332846505024,0.0791227372065423,0.0721591272672233,0.0715330339740876;];

TT=0:60:600;

phi=zeros(11,1);
phi(3)=1;

A=[CC(1:end-1)' phi(2:11) ones(10,1)];
b=CC(2:end)';

ss=A\b;
Csim=CC(1);

for i=2:11
    Csim(i)=ss(1)*Csim(i-1)+ss(2)*phi(i)+ss(3);
end

plot(TT/60,Csim,TT/60,CC,'+')

xlabel('\bfTime [hours]')
ylabel('\bfCreatinine [mmol/L]')
