function [F_true]=cuminc(gamma, alpha,beta1, beta2,beta3,Z,X,delta,V,v,h,n, aipw_beta_hat, nstrt )
% 08/08/2017

n0=length(v);
beta_n=size(Z,2);
K=@(x)0.75*(1-x.^2).*(abs(x)<=1);

z_pd=[0.5,0.5,0.5];
tt=1;

 for j=1:n0  
	hz=(exp(gamma+beta1*z_pd(1)+beta2*z_pd(2)+beta3*z_pd(3))-1)/(gamma+beta1*z_pd(1)+beta2*z_pd(2)+beta3*z_pd(3))*exp(alpha*z_pd(1));
	F_true(j)=exp(alpha*z_pd(1)+(gamma+beta1*z_pd(1)+beta2*z_pd(2)+beta3*z_pd(3))*v(j))/hz*(1-exp(-hz*tt));
end



% Z_pd=z_pd;
% tt=1;

% for j=1:n0

    
      % A0_tx(:,j)=sum((repmat(X,1,n)>=repmat(X',n,1)).*repmat(delta',n,1).* repmat( K((V'-v(j))/h)/h ./sum((repmat(X,1,n)>=repmat(X',n,1))),n,1 ) ,2);
% end
% d_v=diff(v);


% for j=1:n0    
    % eb=0;
    % for k=1:beta_n
        % eb=aipw_beta_hat(j+n0*(k-1))*Z_pd(k)+eb;        
    % end
    % e_A=sum(A0_tx.*d_v(1)*exp(eb),2);
    % F(j)= sum( (X<=tt).*delta.*exp(-e_A).*   K((V-v(j))/h)/h ./   sum((repmat(X,1,n)>=repmat(X',n,1)))' )*exp(eb) ;
% end




%quantile(IMPSTLOG_AUCMB,0.1)
    
 
    