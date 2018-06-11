function [eff_beta_hat,sig_eff_beta,p_values_t1,p_values_t2]=aipw(beta_ini,X,Z,Zt,V,v,delta,sel,beta_Z,beta_Zt,wei,h,n0,n,K,Z_c,method,a,b,aa);
 % 08/08/2017

Z0=Z;
Z0(repmat(sel,1,beta_Z)==0)=NaN;
if ( method=='nonparamet' )
    for k=1:size(Z,2)
        predZ1=ksrmv(Z_c,Z(:,k),repmat(0.2,size(Z_c,2)));
        predZc(:,k)=predZ1.f;
    end
end
if ( method=='parametric' )
    predZ=Z0; 
    Z_1c=[repmat(1,n,1),Z_c];
    for k=1:size(Z,2)
        [theta2,bint,r] = regress(Z0(:,k),Z_1c);
        predZ(:,k)=sum(repmat(theta2',n,1).*Z_1c,2);
    end
end
D=1000;
G=normrnd(0,1,D,n);hat_H=zeros(n,n0);
invA1=zeros(n0,(beta_Z+beta_Zt)^2);

for k=1:beta_Z
    Zc(:,(k-1)*n+1:k*n)=repmat(Z(:,k),1,n);
    predZc(:,(k-1)*n+1:k*n)=repmat(predZ(:,k),1,n);
end
ZZ=[Zc,Zt];
predZZ=[predZc,Zt];
beta_n=beta_Z+beta_Zt;
for k=1:beta_n
    Ze(:,k)=diag(ZZ(:,(k-1)*n+1:k*n));
    predZe(:,k)=diag(predZZ(:,(k-1)*n+1:k*n));
end

for j=1:n0
    beta00=beta_ini(:,j);
    iter=1;Error=1;beta11=beta00;epsilon=0.001;
    while( Error>epsilon & iter < 20)
        eb=0;   predeb=0;  ZZ(repmat(sel,1,beta_Z)==0)=0;
        for k=1:beta_n
            eb=beta00(k)*ZZ(:,((k-1)*n+1):k*n)+eb;
            predeb=beta00(k)*predZZ(:,((k-1)*n+1):k*n)+predeb;
        end

        Jacb=zeros(beta_n,beta_n); AS0=[]; AS1=[]; 
        ZZ(repmat(sel,1,beta_n*n)==0)=0;
        for k=1:beta_n
            AS0(k,:)=mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n));
            AS1(k,:)=mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*ZZ(:,(k-1)*n+1:k*n))+mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n).*predZZ(:,(k-1)*n+1:k*n));
        end

        Ubeta=sum(repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n).*(Ze-(AS1./AS0)')+repmat(K((V-v(j))/h)/h.*delta.*(1-wei),1,beta_n).*(predZe-(AS1./AS0)'));
        for k=1:beta_n
            for l=1:beta_n
                AS2=mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*ZZ(:,(k-1)*n+1:k*n).*ZZ(:,(l-1)*n+1:l*n))+mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n).*predZZ(:,(k-1)*n+1:k*n).*predZZ(:,(l-1)*n+1:l*n));
                Jacb(k,l)=-sum( K((V-v(j))/h)/h.*delta.*wei.* (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))'+K((V-v(j))/h)/h.*delta.*(1-wei).* (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))');
            end
        end
        
        if (rcond(Jacb)<0.00000001)    break;   end
        try  beta11=beta00-pinv(Jacb)*Ubeta';
            Error=sum(abs(beta11-beta00));
            iter=iter+1;
            beta00=beta11;
        catch break;
        end
    end
    eff_beta_hat(1,j:n0:beta_n*n0)=beta11';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%variance estimate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:beta_n
    rr=ksrmv(v',eff_beta_hat((k-1)*length(v)+(1:length(v)))',v(3)-v(1),V);
    eff_beta_u(:,k)=rr.f;
end
eb=0;predeb=0;
for k=1:beta_n
    eb= repmat(eff_beta_u(:,k),1,n).*ZZ(:,((k-1)*n+1):k*n) +eb;
    predeb= repmat(eff_beta_u(:,k),1,n).*predZZ(:,((k-1)*n+1):k*n) +predeb;
end
for k=1:beta_n
    AS0(k,:)=mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n));
    AS1(k,:)=mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*ZZ(:,(k-1)*n+1:k*n))+mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n).*predZZ(:,(k-1)*n+1:k*n));
end
%d lambda
dlambda=1/n./AS0(1,:).*delta'.*sum((repmat(X.*delta,1,n)==repmat(X',n,1)));


for j=1:n0
    Var1=repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n).*(Ze-(AS1./AS0)')+repmat(K((V-v(j))/h)/h.*delta.*(1-wei),1,beta_n).*(predZe-(AS1./AS0)');
    for k=1:beta_n
        Var2(:,k)=-nansum(repmat(K((V-v(j))/h)/h.*wei.*(Ze(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(dlambda,n,1),2)-nansum(repmat(K((V-v(j))/h)/h.*(1-wei).*(predZe(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(dlambda,n,1),2);
    end
    
    for k=1:beta_n
        for l=1:beta_n
            AS2=nanmean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*ZZ(:,(k-1)*n+1:k*n).*ZZ(:,(l-1)*n+1:l*n))+mean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n).*predZZ(:,(k-1)*n+1:k*n).*predZZ(:,(l-1)*n+1:l*n));
            A(k,l)=-nansum( K((V-v(j))/h)/h.*delta.*wei.* (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))'+K((V-v(j))/h)/h.*delta.*(1-wei).* (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))');
        end
    end
    
    sig_eff_beta(1,j:n0:beta_n*n0)=sqrt(diag(pinv(A)'*((Var1+Var2)'*(Var1+Var2))*pinv(A)));
    
    
    
    Var11=repmat(delta.*wei,1,beta_n).*(Ze-(AS1./AS0)')+repmat(delta.*(1-wei),1,beta_n).*(predZe-(AS1./AS0)');
    for k=1:beta_n
        Var22(:,k)=-sum(repmat(wei.*(Ze(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(dlambda,n,1),2)-sum(repmat((1-wei).*(predZe(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(dlambda,n,1),2);
      
    end
    
    invA=pinv(A/n);
    invA1(j,:)=invA(:);
    for ii=1:n
        H_uv=sum(repmat((K((V(ii)-v)/h)/h.*(v>=a).*(v<=v(j)))',1,beta_n^2).*invA1,1)*(v(3)-v(2));
        Varr=(Var11+Var22);
        HH_uv(ii,:)=H_uv;
        hat_H0=reshape(H_uv,beta_n,beta_n)*(Varr(ii,:))';
		for k=1:beta_n
			hat_H(ii,j,k)=hat_H0(k);
		end
    end
    cov_Wb1(j)=1/n*hat_H(:,j)'*hat_H(:,j);
end


%%%%%%%%%%%%test 1%%%%%%%
for k=1:beta_n
    
    eff_beta_hat1=eff_beta_hat';
    eff_beta_hat1(v<=a|v>=b)=0;
    Beta_hat0=[tril(ones(n0))*eff_beta_hat1((k-1)*n0+1:k*n0,:)*(v(3)-v(2))];
    Beta_hat1=sqrt(n)*Beta_hat0;
    Beta_hat=Beta_hat1(1:n0);   %Q^{(1)}
    
    vv=repmat(v,D,1);
    Wb0=G*hat_H(:,:,k)/sqrt(n); Wb=Wb0.*(vv>a&vv<b);
    cov_Wb=[0,cov_Wb1];
    
    
    Beta_hat(Beta_hat==0)=nan;
    Wb(Wb==0)=nan;
    A1_hat=max(abs(Beta_hat));
    A1_tuta=max(abs(Wb),[],2);
    
    A2_hat=nansum(Beta_hat'.^2.*(cov_Wb(2:n0+1)-cov_Wb(1:n0)));
    A2_tuta=nansum(Wb.^2.*repmat((cov_Wb(2:n0+1)-cov_Wb(1:n0)),D,1),2);
    
    M1_hat=min(Beta_hat);
    M1_tuta=min(Wb,[],2);
    
    M2_hat=nansum(Beta_hat'.*(cov_Wb(2:n0+1)-cov_Wb(1:n0)));
    M2_tuta=nansum(Wb.*repmat((cov_Wb(2:n0+1)-cov_Wb(1:n0)),D,1),2);
    
    M3_hat=max(Beta_hat);
    M3_tuta=max(Wb,[],2);
    
    pv_a1(k)=mean(A1_tuta>A1_hat);
    pv_a2(k)=mean(A2_tuta>A2_hat);
    pv_m1(k)=mean(M1_tuta<M1_hat);
    pv_m2(k)=mean(M2_tuta<M2_hat);
    pv_m3(k)=mean(M3_tuta>M3_hat);
    pv_m4(k)=mean(M2_tuta>M2_hat);
    
    
    %%%%%%%%%%%%%%%%%test 2%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    eff_beta_hat1=eff_beta_hat';
    eff_beta_hat1(v<=a|v>=b)=0;
    Beta_hat0=[tril(ones(n0))*eff_beta_hat1((k-1)*n0+1:k*n0,:)*(v(3)-v(2))];
    fir=find(v>a, 1, 'first'); las=find(v<b, 1, 'last');
    Beta_hat1=(Beta_hat0-Beta_hat0(fir))./(v-a)'-(Beta_hat0(las)-Beta_hat0(fir))/(b-a);%a cant be grid point
    Beta_hat2=sqrt(n)*Beta_hat1;
    Beta_hat=Beta_hat2.*(v>aa&v<=b)';   %Q^{(1)}
    Beta_hat(isnan(Beta_hat))=0;
    
    vv=repmat(v,D,1);
    Wb0=G*hat_H(:,:,k)/sqrt(n);
    Wb1=(Wb0-repmat(Wb0(:,fir),1,n0))./repmat((v-a),D,1)-repmat((Wb0(:,las)-Wb0(:,fir))/(b-a),1,n0);
    Wb=Wb1.*(vv>aa&vv<=b);
    Wb(isnan(Wb))=0;
    cov_Wb=[0,cov_Wb1];
    
    
    Beta_hat(Beta_hat==0)=nan;
    Wb(Wb==0)=nan;
    A1_hat=max(abs(Beta_hat));
    A1_tuta=max(abs(Wb),[],2);
    
    A2_hat=nansum(Beta_hat'.^2.*(cov_Wb(2:n0+1)-cov_Wb(1:n0)));
    A2_tuta=nansum(Wb.^2.*repmat((cov_Wb(2:n0+1)-cov_Wb(1:n0)),D,1),2);
    
    M1_hat=min(Beta_hat);
    M1_tuta=min(Wb,[],2);
    
    M2_hat=nansum(Beta_hat'.*(cov_Wb(2:n0+1)-cov_Wb(1:n0)));
    M2_tuta=nansum(Wb.*repmat((cov_Wb(2:n0+1)-cov_Wb(1:n0)),D,1),2);
    
    M3_hat=max(Beta_hat);
    M3_tuta=max(Wb,[],2);
    
    
    pv_a1_t2(k)=mean(A1_tuta>A1_hat);
    pv_a2_t2(k)=mean(A2_tuta>A2_hat);
    pv_m1_t2(k)=mean(M1_tuta<M1_hat);
    pv_m2_t2(k)=mean(M2_tuta<M2_hat);
    pv_m3_t2(k)=mean(M3_tuta>M3_hat);
    pv_m4_t2(k)=mean(M2_tuta>M2_hat);
    
end

p_values_t1=[pv_a1;pv_a2;pv_m1;pv_m2;pv_m3;pv_m4];
p_values_t2=[pv_a1_t2;pv_a2_t2;pv_m1_t2;pv_m2_t2;pv_m3;pv_m4_t2];

