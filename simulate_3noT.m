function [Z,Zt,Aux,X,V,delta ]=simulate_3noT(n,alpha,gamma,beta1,beta2, beta3,strata)
% 08/08/2017

if ( length(gamma)==2 ) gamma1=gamma(1)*strata+gamma(2)*(1-strata);
elseif ( length(gamma)==1 ) gamma1=gamma(1)*ones(n,1);
end

z1=unifrnd(0,1,n,1);
z2=unifrnd(0,1,n,1);
Aux=z1+normrnd(0,0.11,n,1);%auxillary variable  .7

z3=unifrnd(0,1,n,1);

Z=[z1,z2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%generate time and mark%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=exp(-alpha*z1).*(gamma1+beta1*z1+beta2*z2+beta3*z3)./(exp(gamma1+beta1*z1+beta2*z2+beta3*z3)-1);
failtime=-log(unifrnd(0,1,n,1)).*mu;
censortime=min(exprnd(4,n,1),2);

delta=(failtime<=censortime);
X=min(failtime,censortime);

%V=log(unifrnd(0,1,n,1).*(exp(gamma+beta1*z1+beta2*z2+beta3*z3)-1)+1)./(gamma+beta1*z1+beta2*z2+beta3*z3);%mark
V=log(unifrnd(0,1,n,1).*(exp(gamma1+beta1*z1+beta2*z2+beta3*z3)-1)+1)./(gamma1+beta1*z1+beta2*z2+beta3*z3);%mark

z3t=repmat(z3,1,n);
Zt=[z3t];



