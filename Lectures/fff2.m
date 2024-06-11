function y=fff2(a,t,F)
x=5*(exp(-a(1))+exp(-a(2)))+a(2)*cos(a(:,1).*a(2).^2*t)+a(2)^2*sin(a(1)^2/(a(2)+.1)*t);
y=norm(x-F);
