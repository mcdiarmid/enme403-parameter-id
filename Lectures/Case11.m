clc
beep off
clear all
close all
format short g
set(0,'units','centimeters','defaultfigureposition',[60 60 420,330])

figure(1101)
psi2=zeros(121,101); % it's always a good idea to make these non-square as it becomes easy to plt the resulting matrix agains the wrong axis!

for a1=0:0.01:1.2
    for a2=0:.01:1
        %fill the matrix with the objective value
        psi2(round(a1*100+1),round(a2*100+1))=1.5*a1^2-2.5*a1+2*a1*a2-2.1*a2+a2^2+1.1825;
    end
end

contour(0:0.01:1,0:0.01:1.2,psi2);hold all
plot(0.65,0.4,'+r','linewidth',2,'markersize',10)
xlabel('\bf\it\theta_2');
ylabel('\bf\it\theta_1');

figure(1102)

set(gcf,'position',[60 60 680 660]) %setup some prettier axii than subplot
AA1=axes('position',[0.08 0.57 0.4 0.42]);hold all
xlabel('\bf\it\theta_2');ylabel('\bf\it\theta_1')
AA2=axes('position',[0.58 0.57 0.4 0.42]);hold all
xlabel('\bfiterations');ylabel('\bf\it\theta_i');xlim([0 100])
AA3=axes('position',[0.08 0.07 0.4 0.42]);hold all
xlabel('\bfiterations');ylabel('\bf\it\psi');xlim([0 100])

I=1;
a12=[0  0];
psi=zeros(1,100);
% the model residual at a particular parameter set
psi(1)=1.5*a12(I,1)^2-2.5*a12(I,1)+2*a12(I,1)*a12(I,2)-2.1*a12(I,2)+a12(I,2)^2+1.1825;
% the model Jacobian at a particular parameter set
jj12=[3*a12(1)+2*a12(2)-2.5 2*a12(1)+2*a12(2)-2.1];
contour(AA1,0:0.01:1,0:0.01:1.2,psi2)

plot(AA1,0.65,0.4,'+r','linewidth',2,'markersize',10)

for I=2:100
    a12(I,:)=a12(I-1,:)-jj12(I-1,:)*0.1;        % go downhill
    jj12(I,:)=[3*a12(I,1)+2*a12(I,2)-2.5 2*a12(I,1)+2*a12(I,2)-2.1];
    psi(I)=1.5*a12(I,1)^2-2.5*a12(I,1)+2*a12(I,1)*a12(I,2)-2.1*a12(I,2)+a12(I,2)^2+1.1825;
    
    plot(AA1,a12(:,2),a12(:,1),'.-k')
    plot(AA2,1:I,a12(1:I,1),'.-b');
    plot(AA2,1:I,a12(1:I,2),'.-g');
    semilogy(AA3,1:I,psi(1:I),'.-k');           % log plot!
    drawnow;pause(.25*norm(jj12(I,:)))          % slow down according to the step length!
%     if I==3
%         dsfg jkldsfg hdjkfg dfjk              % better than 'break' because it ALWAYS
%                                               % stops the iterations
%     end
end

legend(AA2,'\bf\it\theta_1','\bf\it\theta_2')             % legends take a long time to plot, strangely...

%% MonteCarlo
figure(1103)
set(gcf,'position',[60 60 680 660]) %setup some prettier axii than subplot
AA1=axes('position',[0.08 0.57 0.4 0.42]);hold all
xlabel('\bf\it\theta_2');ylabel('\bf\it\theta_1')
AA2=axes('position',[0.58 0.57 0.4 0.42]);hold all
xlabel('\bfiterations');ylabel('\bf\it\theta_i');xlim([0 100])
AA3=axes('position',[0.08 0.07 0.4 0.42]);hold all
xlabel('\bfiterations');ylabel('\bfrate (%)');
contour(AA1,0:0.01:1,0:0.01:1.2,psi2)
for III=1:100
    % our I.C.s
    a12=[0.6+0.6*sin(2*pi*III/100) 0.5+0.5*cos(2*pi*III/100)];
    jj12=[3*a12(1)+2*a12(2)-2.5 2*a12(1)+2*a12(2)-2.1];
    
    for I=2:200
        a12(I,:)=a12(I-1,:)-jj12(I-1,:)*0.1;
        jj12(I,:)=[3*a12(I,1)+2*a12(I,2)-2.5 2*a12(I,1)+2*a12(I,2)-2.1];
        psi(I)=1.5*a12(I,1)^2-2.5*a12(I,1)+2*a12(I,1)*a12(I,2)-2.1*a12(I,2)+a12(I,2)^2+1.1825;
        
        plot(AA1,a12(:,2),a12(:,1),'.-k')
        plot(AA2,1:I,a12(1:I,1),'-b');
        plot(AA2,1:I,a12(1:I,2),'-g');
    end
    %find how many iterations were required for convergence
    conv(III)=min([find((abs(a12(:,1)-0.4)<0.004).*(abs(a12(:,2)-0.65)<0.0065)) ;999]);
    
end
plot(AA1,0.65,0.4,'+r','linewidth',2,'markersize',10)
plot(AA3,sort(conv),1:100,'k');grid on

% saveas(gcf,'figggg.fig')
%%
figure(1104);
set(gcf,'position',[60 60 420 660]) %setup some prettier axii than subplot
AA1=axes('position',[0.08 0.57 0.88 0.42]);hold all
xlabel('\bf\it\theta_2');ylabel('\bf\it\theta_1');axis([0 1 0 1.2])
AA2=axes('position',[0.08 0.07 0.88 0.42]);hold all
xlabel('\bfiterations');ylabel('\bfincidence');xlim([0 200])

hold all
contour(AA1,0:0.01:1,0:0.01:1.2,psi2)


for alpha=[.05 .1 .2 .3 .4 .5]
    
    I=1;
    a12D=[0 0];
    dela12=[0.0004 0.00065];
    % I got bored of writing out function so made a function to evaluate psi C11fcn.
    jj12D=[(C11fcn(a12D(I,1)+dela12(1),a12D(I,2))-C11fcn(a12D(I,1),a12D(I,2)))/dela12(1) (C11fcn(a12D(I,1),a12D(I,2)+dela12(2))-C11fcn(a12D(I,1),a12D(I,2)))/dela12(2)];
    
    for I=2:100
        a12D(I,:)=a12D(I-1,:)-jj12D(I-1,:)*alpha;
        jj12D(I,:)=[(C11fcn(a12D(I,1)+dela12(1),a12D(I,2))-C11fcn(a12D(I,1),a12D(I,2)))/dela12(1) (C11fcn(a12D(I,1),a12D(I,2)+dela12(2))-C11fcn(a12D(I,1),a12D(I,2)))/dela12(2)];
    end
    plot(AA1,a12D(:,2),a12D(:,1),'.-');hold all
    
    
end
legend(AA1,{'\it\psi' '.05' '.1' '.2' '.3' '.4' '.5'})

%% MonteCarlo


xx=0;
for alpha=[.05 .1 .2 .3 .4 .5]
    xx=xx+1;
    for III=1:100
        a12=[0.6+0.6*sin(2*pi*III/100) 0.5+0.5*cos(2*pi*III/100)];
        jj12=[3*a12(1)+2*a12(2)-2.5 2*a12(1)+2*a12(2)-2.1];
                
        for I=2:200
            a12(I,:)=a12(I-1,:)-jj12(I-1,:)*alpha;
            jj12(I,:)=[3*a12(I,1)+2*a12(I,2)-2.5 2*a12(I,1)+2*a12(I,2)-2.1];
        end
        %find how many iterations were required for convergence
        conv(III)=min([find((abs(a12(:,1)-0.4)<0.004).*(abs(a12(:,2)-0.65)<0.0065)) ;999]);
        
    end
    plot(AA2,sort(conv),1:100,'-');hold all
end

xlim([0 200])


%% Discrete
figure(1105);hold all
                        
I=1;
a12=[0  0];
dela=[0.0004 0.00065];  % some randomly chosen perterbations
psi=zeros(1,100);
psi(1)=C11fcn(a12(1),a12(2));

jj12=[(C11fcn(a12(I,1)+dela(1),a12(I,2))-C11fcn(a12(I,1),a12(I,2)))/dela(1) (C11fcn(a12(I,1),a12(I,2)+dela(2))-C11fcn(a12(I,1),a12(I,2)))/dela(2) ];
contour(0:0.01:1,0:0.01:1.2,psi2);hold all

plot(0.65,0.4,'+r','linewidth',2,'markersize',10)

for I=2:100 % all else is the same
    a12(I,:)=a12(I-1,:)-jj12(I-1,:)*0.1;
    jj12(I,:)=[(C11fcn(a12(I,1)+dela(1),a12(I,2))-C11fcn(a12(I,1),a12(I,2)))/dela(1) (C11fcn(a12(I,1),a12(I,2)+dela(2))-C11fcn(a12(I,1),a12(I,2)))/dela(2) ];
    plot(a12(:,2),a12(:,1),'.-k')
end
xlabel('\bf\it\theta_2');ylabel('\bf\it\theta_1');axis([ 0 1 0 1.2])
