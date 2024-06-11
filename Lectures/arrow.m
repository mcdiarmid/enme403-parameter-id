function angle=arrow2(X,Y,a,AR)
% pdd 13-5-11
% X=X1-X2
% Y=Y1-Y2
% a=color i.e. 'k'
% AR aspect ratio of figure axis - if it doens work try AR=1/AR
X
Y


XY=[X(:) Y(:)]'
X=XY(:,1)
Y=XY(:,2)

plot(X,Y,a);hold all
arrowhead=[0 0;5/6*pi .2;7/6*pi 0.2];
arrowhead(:,2)=arrowhead(:,2)*((X(2)-X(1))^2+(Y(2)-Y(1))^2)^0.5;

angle=atan((X(2)-X(1))/(Y(2)-Y(1)));
if Y(2)<Y(1)
    angle=angle+pi;
end
    
arrowhead(:,1)=arrowhead(:,1)+angle;
arr=([sin(arrowhead(1,1))*arrowhead(1,2) sin(arrowhead(2,1))*arrowhead(2,2) sin(arrowhead(3,1))*arrowhead(3,2);cos(arrowhead(1,1))*arrowhead(1,2) cos(arrowhead(2,1))*arrowhead(2,2) cos(arrowhead(3,1))*arrowhead(3,2)]');

mnn=(arr(2,1)+arr(3,1))/2;
dif=(arr(2,1)-arr(3,1))/2;

arr(2,1)=mnn+dif*AR;
arr(3,1)=mnn-dif*AR;

mnn=(arr(2,2)+arr(3,2))/2;
dif=(arr(2,2)-arr(3,2))/2;

arr(2,2)=mnn+dif/AR;
arr(3,2)=mnn-dif/AR;

arr(:,1)=arr(:,1)+X(2);
arr(:,2)=arr(:,2)+Y(2);

fill(arr(:,1),arr(:,2),a)

