function CCt=CCtsim(k,U,t,UX,V,C0)
CCt=exp(-k*t).*(C0+cumtrapz(t,exp(k*t).*(UX+U)/V));