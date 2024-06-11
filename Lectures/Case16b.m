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

figure(1611)
contour(aa,aa,psi,50);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')

orgs=5*rand(100,3);
plot(orgs(:,3),orgs(:,2),'.k')

figure(1614)
contour(aa,aa,psi,50);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')

plot(orgs(:,3),orgs(:,2),'.k')

figure (1612)
contour(aa,aa,psi,50);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')

plot(orgs(:,3),orgs(:,2),'.k')


for generations=1:10
    
    
    for I=1:length(orgs)
        orgs(I,1)=norm(fff(orgs(I,2),orgs(I,3),t)-F);
        
    end
    orgs=sortrows(orgs);
    maxerr(generations,:)=[orgs(1,:) orgs(end,1)];
    
    % plot3(orgs(:,3),orgs(:,2),orgs(:,1),'.k');hold all;drawnow
    plot(orgs(:,3),orgs(:,2),'.k');hold all;drawnow;pause(0.5)
    
    orgs(:,2:3)=orgs(ceil(((1:100)/14).^2),2:3);
    
    orgs(2:end,2:3)=orgs(2:end,2:3)+0.01*[orgs(2:end,1) orgs(2:end,1)].*randn(99,2);
    % orgs(2:end,2:3)=orgs(2:end,2:3)+0.2*randn(99,2);
%     maxerr(generations,:)=[orgs(1,:) orgs(end,1)];
    % orgs(2:end,2:3)=max(0,min(5,(orgs(2:end,2:3))));
    
end
figure(1614)
plot(maxerr(:,3),maxerr(:,2),'-r','linewidth',3)


figure(1613)
subplot(2,1,1);plot(maxerr(:,2:3));xlabel('\bfgenerations');ylabel('\bf\itx')
subplot(2,1,2);plot(maxerr(:,1));xlabel('\bfgenerations');ylabel('\bf\it\psi')










