function CT=fCCt(k,Un,t,Ux,C0,V,T)
% just gives data given a particular parameter set
C=exp(-k*t).*(C0+cumtrapz(t,exp(k*t).*(Ux+Un)/V));
CT=C(T);
