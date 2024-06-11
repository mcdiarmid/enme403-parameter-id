clc
beep off
clear all
close all
format short g
warning off

set(0,'units','centimeters','defaultfigureposition',[60 60 420 330])

t=0:.1:12;

% I needed 4 roots so this is how I set up the model:
x=(t-1).*(t-5).*(t-6).*(t-7)+65;
figure(1501);
plot(t,x);
xlabel('\bf\itt');ylabel('\bf\itx');ylim([0 300])
a=polyfit(t,x,4);

                        
for IC=0:.1:12          % cycle through initial conditions
    v=zeros(50,1);
    v(1:2)=[999 IC];
    I=2;                            
    while abs(v(I)-v(I-1))>0.001    % keep going till you stop
        I=I+1;
        v(I)=v(I-1)-0.01*(4*a(1)*v(I-1)^3+3*a(2)*v(I-1)^2+2*a(3)*v(I-1)+a(4));
%         figure(15011)
%         plot(t,x);hold all
%         plot(v(I),polyval(a,v(I)),'+');drawnow;%pause(0.1)
%         hold off
       v(I)=min(12,max(0,v(I)));
    end
    v(I:end)=v(I);
    figure(1502)
    plot(v,'.-k');hold all
    ICo(round(10*IC+1))=v(I);
    axis([0 40 0 12])
end
xlabel('\bfiterations');ylabel('\bfargmin\it_tx(t)')
figure(1503)
plot(0:.1:12,ICo,'.')
xlabel('\bfI.C.');ylabel('\bfargmin\it_tx(t)')

figure(1504)
plot(t,x);hold all
xlabel('\bf\itt');ylabel('\bf\itx');ylim([0 300])

x=(0:10)';
for I=1:10
    psi=(x(:,1)-1).*(x(:,1)-5).*(x(:,1)-6).*(x(:,1)-7)+65; %evaluate the fitness
    x(:,2)=psi;
    plot(x(:,1),x(:,2),'.k');drawnow
    x=sortrows(x,2);                                        %order according to fitness
    %save the best and get some variance about that and other good ones
    x(:,1)=[x(1,1);x(1,1)+0.5*randn(4,1);x(2,1)+0.5*randn(3,1);x(3,1)+0.5*randn(2,1);x(4,1)+0.5*randn(1,1)];
    
    
end

x=x(1,1)



















