function [Z,Zt,Aux,X,V,delta,strata ]=simulate2(n,alpha,gamma,beta1,beta2, beta3,strata)
% 08/08/2017

if ( length(gamma)==2 ) gamma1=gamma(1)*strata+gamma(2)*(1-strata);
elseif ( length(gamma)==1 ) gamma1=gamma(1)*ones(n,1);
end

z1=unifrnd(0,1,n,1)/2 ;
u2=unifrnd(0,1,n,1);
z2=(-0.1*z1+u2)/2;
 
Aux=z1+normrnd(0,0.15,n,1);  %auxillary variable  .7

z3=unifrnd(0,1,n,1)/4;

Z=[z1,z2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%generate time and mark%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bz=exp(alpha*z1).*(exp(gamma1+beta1*z1+beta2*z2)-1)./(gamma1+beta1*z1+beta2*z2);
failtime=(-(beta3*z3+1).*log(unifrnd(0,1,n,1))./bz).^(1./(beta3*z3+1));

censortime=min(exprnd(4,n,1),2);
delta=(failtime<=censortime);
X=min(failtime,censortime);

V=(log(unifrnd(0,1,n,1).*(exp(gamma1+beta1*z1+beta2*z2)-1)+1)./(gamma1+beta1*z1+beta2*z2)).^(1/2); %marks

z3t=repmat(z3,1,n).*log(repmat(failtime',n,1));
Zt=[z3t];