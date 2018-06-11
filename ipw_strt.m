function [beta_hat,sig_beta]=ipw_strt(beta_ini,X,Z,Zt,V,v,delta,sel,beta_Z,beta_Zt,wei,h,n0,n,K,strata)
 % 08/08/2017
nstrt=length(unique([strata]));
ZZ=Z;

for k=1:beta_Z
    Z_all1(:,(k-1)*n+1:k*n)=repmat(Z(:,k),1,n);
end
Z_all=[Z_all1,Zt];
beta_n=beta_Z+beta_Zt; Ze=[];
for k=1:beta_n
    Ze(:,k)=diag(Z_all(:,(k-1)*n+1:k*n));
end

for j=1:n0
    beta00=beta_ini(:,j);
    iter=1;Error=1;beta11=beta00;epsilon=0.001;
    while( Error>epsilon & iter < 50)
        eb=0;
        for k=1:beta_n
            eb=beta00(k)*Z_all(:,((k-1)*n+1):k*n)+eb;
        end
        for kk=1:nstrt
            S0kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))./sum(repmat(strata==kk-1,1,n));
            S0kk1(kk,:)=S0kk.*(strata==kk-1)';
        end
        S0=sum(S0kk1,1);
        
        Jacb=zeros(beta_n,beta_n);
        for k=1:beta_n
            for kk=1:nstrt
                S1kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*Z_all(:,((k-1)*n+1):k*n))./sum(repmat(strata==kk-1,1,n));
                S1kk1(kk,:)=S1kk.*(strata==kk-1)';
            end
            S1(k,:)=sum(S1kk1,1);
        end
        Ubeta=sum(repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n).*(Ze-(S1./(repmat(S0,beta_n,1)))'));
        for k=1:beta_n
            for l=1:beta_n
                for kk=1:nstrt
                    S2kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*Z_all(:,(k-1)*n+1:k*n).*Z_all(:,(l-1)*n+1:l*n))./sum(repmat(strata==kk-1,1,n));
                    S2kk1(kk,:)=S2kk.*(strata==kk-1)';
                end
                S2=sum(S2kk1,1);
                Jacb(k,l)=-sum( K((V-v(j))/h)/h.*delta.*wei.* (S2./S0-(S1(k,:)./S0.*S1(l,:)./S0))');
            end
        end
        
        if (rcond(Jacb)<0.00000001)    break;   end
        try  beta11=beta00-pinv(Jacb)*Ubeta';
            Error=sum(abs(beta11-beta00)) ;
            iter=iter+1;
            beta00=beta11;
        catch break;
        end
    end
    iter;
    beta_hat(1,j:n0:beta_n*n0)=beta11';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%variance estimate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eb=0;
    for k=1:size(Z,2)
        eb=beta11(k)*Z_all(:,((k-1)*n+1):k*n)+eb;
    end
    for kk=1:nstrt
        S0kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))./sum(repmat(strata==kk-1,1,n));
        S0kk1(kk,:)=S0kk.*(strata==kk-1)';
    end
    S0=sum(S0kk1,1);
    %d lambda
   
    Jacb=zeros(beta_n,beta_n);
    for k=1:beta_n
        for kk=1:nstrt
            S1kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*Z_all(:,((k-1)*n+1):k*n))./sum(repmat(strata==kk-1,1,n));
            S1kk1(kk,:)=S1kk.*(strata==kk-1)';
        end
        S1(k,:)=sum(S1kk1,1);
    end
    Var1=repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n).*(Ze-(S1./(repmat(S0,beta_n,1)))');
    for k=1:beta_n
        for kk=1:nstrt
            S0kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*repmat(strata==kk-1,1,n).*repmat(strata==kk-1,1,n)')./sum(repmat(strata==kk-1,1,n));
            dlambdakk=1./S0kk.*sum((repmat(X.*(strata==kk-1).*delta,1,n)==repmat(X',n,1)).*repmat(strata==kk-1,1,n)')./sum(repmat(strata==kk-1,1,n));
            dlambdakk(isnan(dlambdakk))=0;
            Var2a=-sum(repmat(K((V-v(j))/h)/h.*wei.*(Ze(:,k)-(S1(k,:)./S0)'),1,n).*(repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*repmat(strata==kk-1,1,n).*exp(eb).*repmat(dlambdakk,n,1),2);
            Var2kk(kk,:)=Var2a'.*(strata==kk-1)';
        end
        Var2(:,k)=sum(Var2kk,1)';
    end
    
    for k=1:beta_n
        for l=1:beta_n
            for kk=1:nstrt
                S2kk=sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*Z_all(:,(k-1)*n+1:k*n).*Z_all(:,(l-1)*n+1:l*n).*repmat(strata==kk-1,1,n))./sum(repmat(strata==kk-1,1,n));
                S2kk1(kk,:)=S2kk.*(strata==kk-1)';
            end
            S2=sum(S2kk1,1);
            A(k,l)=-sum( K((V-v(j))/h)/h.*delta.*wei.* (S2./S0-(S1(k,:)./S0.*S1(l,:)./S0))');
        end
    end
    sig_beta(1,j:n0:beta_n*n0)=sqrt(diag(pinv(A)'*((Var1+Var2)'*(Var1+Var2))*pinv(A)));
end