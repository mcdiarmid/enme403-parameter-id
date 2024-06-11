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

figure(1601)
surf(aa,aa,psi,'linestyle','none');hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')
zlabel('\bf\it\psi')

% plot a bunch of starting conditions
x=meshgrid(0:.25:5,0:.25:5);
plot3(x,x',interp2(aa,aa,psi,x,x'),'.k','markerfacecolor','k')
% plot(x(1,:),x(2,:),'.k')
view([-20 20])          % choose some sweet angle

figure(1602)
contour(aa,aa,psi,50);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')
zlabel('\bf\it\psi')
plot(x,x','.k','markerfacecolor','k')


figure(1603)
contour(aa,aa,psi,100);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')
plot(2.5,2.5,'+r')

load matlab.mat

imagesc([0 5],[5 0],(mmap2(end:-1:1,:))*64/24);hold all
colormap(lines)

contour(aa,aa,psi/20000,100,'linecolor','w');hold all

plot(b1,a1,'+r','linewidth',2,'markersize',10)
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')
x=zeros(250,2);

fcs=zeros(100,6);
lambda=5;
ind=0;
for ICx=0:.25:5            % try 100 initial conditions
    for ICy=0:.25:5
        ind=ind+1;
        x(1,:)=[ICx,ICy];
        
        for its=2:500
            [J,psii]=dfff(x(its-1,1),x(its-1,2),t,F);
            dx=(J'*psii')';
            if its<50
                x(its,:)=x(its-1,:)-0.1*dx/max([1, abs(dx)]);
                
            elseif its<100
                x(its,:)=x(its-1,:)-0.05*dx/max([1, abs(dx)]);
                
            elseif its<150
                x(its,:)=x(its-1,:)-0.025*dx/max([1, abs(dx)]);
                
            else
                x(its,:)=x(its-1,:)-0.01*dx/max([1, abs(dx)]);
            end
            
        end
        plot(x(:,2),x(:,1),'-k');hold all
        plot(x(1,2),x(1,1),'.k','markerfacecolor','k');hold all
        plot(x(end,2),x(end,1),'ok','markerfacecolor','w');hold all
        fcs(ind,:)=[x(end,:) norm(psii) 1 x(1,:)];
        solmap(ICx*4+1,ICy*4+1,1:4)=[x(end,:) norm(psii) 1];
    end
end

for ii=1:21
    for jj=1:21
        solmap(ii,jj,4)=0.5*round(2*solmap(ii,jj,1))*5+0.5*round(2*solmap(ii,jj,2))+50;
    end
end

mmap=squeeze(solmap(:,:,4));
mmap2=mmap;
mmap3=mmap(:)
no=1;
for I=50:.5:100
    if find(mmap==I)
        mmap2(find(mmap==I))=no;
        [xx,yy]=find(mmap==I,1,'first');
        fcss(no,:)=[length(find(mmap==I)) solmap(xx,yy,1) solmap(xx,yy,2) round(solmap(xx,yy,3)*10)/10]
        no=no+1;
        ress(no,:)=[xx yy ];
        
    end
end
axis([-.1 5.1 -.1 5.1])

figure(1604)

contour(aa,aa,psi,50);hold all
xlabel('\bf\itx_2')
ylabel('\bf\itx_1')
zlabel('\bf\it\psi')


for I=1:length(fcss)
    
    text(fcss(I,3),fcss(I,2),['\it\bf\psi: ' num2str(fcss(I,4)) '\rm\bf n:' num2str(fcss(I,1))],'backgroundcolor','w')
    plot(fcss(I,3),fcss(I,2),'+k')
end

figure(1605)
psis=round(solmap(:,:,3)/2)*2;
plot(sort(psis(:)),'.');
xlabel('\bfrate')
ylabel('\bf\it\psi')

figure(1606)
fcss=sortrows(fcss,4);

plot(t,F,'k+');hold all
for I=1:5
    plot(t,fff(fcss(I,2),fcss(I,3),t))
end
xlabel('\bf\ittime')
ylabel('\bf\itsome property!')
legend('data','\psi1','\psi2','\psi3','\psi4','\psi5')
% imagesc([0 5],[0 5],mmap2)
% colormap(lines)













