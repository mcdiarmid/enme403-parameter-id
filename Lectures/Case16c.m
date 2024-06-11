clc
beep off
clear all
close all
format short g
warning off

set(0,'units','centimeters','defaultfigureposition',[60 60 680 660])

t = 0:.02:1;

U=0*t;
U(100:1000:end)=1;

%pick some parameter values:
a1=2.5;
b1=2.5;
F=fff(a1,b1,t); %call some function to get 'data', F

xx=0;
aa=linspace(0,5,101);
for a=aa
    xx=xx+1;yy=0;
    for b=aa
        yy=yy+1;
        psi(xx,yy)=norm(fff(a,b,t)-F);
    end
end

figure(1621)
contour(aa,aa,psi,50);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')

x=5*rand(2,3);
plot(x(:,2),x(:,1),'.k')
x(1,3)=norm(fff(x(1,1),x(1,2),t)-F);
plot(x(1,2),x(1,1),'+r');drawnow
maxerr=zeros(1000,3);
for generations=1:1000;
    
    x(2,:)=x(1,:)+0.05*x(1,3)*randn(1,3);
    x(2,1:2)=max(0,min(5,x(2,1:2)));
    plot(x(2,2),x(2,1),'.k');drawnow
    x(2,3)=norm(fff(x(2,1),x(2,2),t)-F);
    
    if x(2,3)<x(1,3)
        x(1,:)=x(2,:);
        plot(x(1,2),x(1,1),'+r');drawnow
        
    end
    maxerr(generations,1:3)=x(1,:);
end

plot(maxerr(:,2),maxerr(:,1),'r','linewidth',2)

figure(1622)
subplot(2,1,1);plot(maxerr(:,1:2));xlabel('\bfgenerations');ylabel('\bf\itx')
subplot(2,1,2);plot(maxerr(:,3));xlabel('\bfgenerations');ylabel('\bf\it\psi')

