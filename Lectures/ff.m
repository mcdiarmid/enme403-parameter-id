function Tt=ff(x,JJ,Tamb,t)
global V cP
Tt=exp(-x(1)*t).*(Tamb+cumtrapz(t,exp(x(1)*t).*(x(1)*Tamb+1000/V/cP/t(2)*(JJ*x(2:4))')));