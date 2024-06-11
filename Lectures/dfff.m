function [J,psi]=dfff(a,b,t,x)
fab=fff(a,b,t);
fda=fff(a+0.001,b,t);
fdb=fff(a,b+0.001,t);
psi=fab-x;
% plot([x;fab;fda;fdb]');


J=[fda-fab;fdb-fab]'/0.001;


