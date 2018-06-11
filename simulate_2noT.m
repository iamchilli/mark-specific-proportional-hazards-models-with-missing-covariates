function [Z,Aux,X,V,delta ]=simulate_2noT(n,alpha,gamma,beta1,beta2,aux_corr)

%strata=binornd(1,0.5,n,1);
z1=unifrnd(0,1,n,1);
z2=binornd(1,0.5,n,1);

 
if (aux_corr==0)  Aux=normrnd(0,0.3,n,1); end  %auxillary variable 0
if (aux_corr==0.7)  Aux=z1+normrnd(0,0.3,n,1); end %auxillary variable  .7
if (aux_corr==0.93)  Aux=z1+normrnd(0,0.11,n,1); end %auxillary variable .93
 
Z=[z1,z2];
beta_n=size(Z,2);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%generate time and mark%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu=exp(-alpha*z1).*(gamma+beta1*z1+beta2*z2+beta3*z3)./(exp(gamma+beta1*z1+beta2*z2+beta3*z3)-1);
mu=exp(-alpha*z1).*(gamma+beta1*z1+beta2*z2)./(exp(gamma+beta1*z1+beta2*z2)-1);
failtime=-log(unifrnd(0,1,n,1)).*mu;
censortime=min(exprnd(4,n,1),2);

delta=(failtime<=censortime);
X=min(failtime,censortime);

%V=log(unifrnd(0,1,n,1).*(exp(gamma+beta1*z1+beta2*z2+beta3*z3)-1)+1)./(gamma+beta1*z1+beta2*z2+beta3*z3);%mark
V=log(unifrnd(0,1,n,1).*(exp(gamma+beta1*z1+beta2*z2)-1)+1)./(gamma+beta1*z1+beta2*z2);%mark
 
