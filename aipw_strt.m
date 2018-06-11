function [eff_beta_hat,sig_eff_beta,p_values_t1,p_values_t2,trace1_sim,trace1_real,trace2_sim,trace2_real, F]=aipw_strt(beta_ini,X,Z,Zt,V,v,delta,sel,beta_Z,beta_Zt,wei,h,n0,n,K,Z_c,method,a,b,aa,strata, Z_pd, tt);
% 08/08/2017
nstrt=length(unique([strata]));
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
G=normrnd(0,1,D,n);
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
hat_H=zeros([n,n0,beta_n]);

%%%%%%%%%%%%%%%% By grid points of v%%%%%%%%%%%%%%%%%%

for j=1:n0
    
    %%%%%%%%%%%%%%%% Estimation %%%%%%%%%%%%%%%%%%
    
    if (j==1) beta00=beta_ini(:,j);
    else beta00=beta11;
    end
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
            for kk=1:nstrt
                S0kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n)))./sum(repmat(strata==kk-1,1,n));
                S1kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*ZZ(:,(k-1)*n+1:k*n))+nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n).*predZZ(:,(k-1)*n+1:k*n)))./sum(repmat(strata==kk-1,1,n));
                S0kk1(kk,:)=S0kk.*(strata==kk-1)';
                S1kk1(kk,:)=S1kk.*(strata==kk-1)';
            end
            AS0(k,:)=sum(S0kk1,1);
            AS1(k,:)=sum(S1kk1,1);
        end
        
        Ubeta=nansum(repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n).*(Ze-(AS1./AS0)')+repmat(K((V-v(j))/h)/h.*delta.*(1-wei),1,beta_n).*(predZe-(AS1./AS0)'));
        for k=1:beta_n
            for l=1:beta_n
                for kk=1:nstrt
                    S2kk=(sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n).*ZZ(:,(k-1)*n+1:k*n).*ZZ(:,(l-1)*n+1:l*n))+sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n).*predZZ(:,(k-1)*n+1:k*n).*predZZ(:,(l-1)*n+1:l*n)))./sum(repmat(strata==kk-1,1,n));
                    S2kk1(kk,:)=S2kk.*(strata==kk-1)';
                end
                AS2=sum(S2kk1,1);
                Jacb(k,l)=-nansum( K((V-v(j))/h)/h.*delta.*wei.* (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))'+K((V-v(j))/h)/h.*delta.*(1-wei).* (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))');
                
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
    
    %%%%%%%%%%%%%%%% Variance estimate%%%%%%%%%%%%%%%%%%
    
    Var1=repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n).*(Ze-(AS1./AS0)')+repmat(K((V-v(j))/h)/h.*delta.*(1-wei),1,beta_n).*(predZe-(AS1./AS0)');
    Var1(isnan(Var1))=0;
    for k=1:beta_n
        for kk=1:nstrt
            S0kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n)))./nansum(repmat(strata==kk-1,1,n));
            dlambdakk=1./S0kk.*delta'.*sum((repmat(X.*delta.*(strata==kk-1),1,n)==repmat(X',n,1)).*repmat(strata==kk-1,1,n)')./sum(repmat(strata==kk-1,1,n));
            Var2a=-nansum(repmat(K((V-v(j))/h)/h.*wei.*(Ze(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(dlambdakk,n,1),2)-nansum(repmat(K((V-v(j))/h)/h.*(1-wei).*(predZe(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(dlambdakk,n,1),2);
            Var2kk(kk,:)=Var2a'.*(strata==kk-1)' ;            
        end
        Var2(:,k)=sum(Var2kk,1)';
    end
    
    A=Jacb;    
    sig_eff_beta(1,j:n0:beta_n*n0)=sqrt(diag(pinv(A)'*((Var1+Var2)'*(Var1+Var2))*pinv(A)));    
    Var11=repmat(delta.*wei,1,beta_n).*(Ze-(AS1./AS0)')+repmat(delta.*(1-wei),1,beta_n).*(predZe-(AS1./AS0)');    
    for k=1:beta_n
        for kk=1:nstrt
            S0kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n)))./sum(repmat(strata==kk-1,1,n));
            dlambdakk=1./S0kk.*delta'.*nansum((repmat(X.*(strata==kk-1).*delta,1,n)==repmat(X',n,1)).*repmat(strata==kk-1,1,n)')./nansum(repmat(strata==kk-1,1,n));
            Var2a=-nansum(repmat( wei.*(Ze(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(dlambdakk,n,1),2)-nansum(repmat( (1-wei).*(predZe(:,k)-(AS1(k,:)./AS0(k,:))'),1,n).*(repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(dlambdakk,n,1),2);
            Var2kk(kk,:)=Var2a'.*(strata==kk-1)';
        end
        Var22(:,k)=nansum(Var2kk,1)';
    end    
    
    invA=pinv(A/n);
    invA1(j,:)=invA(:);
    for ii=1:n
        H_uv=nansum(repmat((K((V(ii)-v)/h)/h.*(v>=a).*(v<=v(j)))',1,beta_n^2).*invA1,1)*(v(3)-v(2));
        Varr=(Var11+Var22);
        HH_uv(ii,:)=H_uv;
        hat_H0=reshape(H_uv,beta_n,beta_n)*(Varr(ii,:))';
        for k=1:beta_n
            hat_H(ii,j,k)=hat_H0(k);
        end
    end
    for k=1:beta_n
        cov_Wb1(j,k)=1/n*hat_H(:,j,k)'*hat_H(:,j,k);
    end
        
end

%%%%%%%%%%%%%%%% Test 1 %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Test 1 %%%%%%%%%%%%%%%%%%

for k=1:beta_n    
    eff_beta_hat1=eff_beta_hat((k-1)*n0+1:k*n0)';
    eff_beta_hat1(v<a|v>b)=0;
    Beta_hat0=[tril(ones(n0))*eff_beta_hat1*(v(3)-v(2))];
    Beta_hat1=sqrt(n)*Beta_hat0;
    Beta_hat=Beta_hat1(1:n0);   %Q^{(1)}
    
    vv=repmat(v,D,1);
    Wb0=G*hat_H(:,:,k)/sqrt(n);
    Wb=Wb0.*(vv>=a&vv<=b);
    Wb(isnan(Wb))=0;
    cov_Wb=[0,cov_Wb1(:,k)'];
        
    Beta_hat(Beta_hat==0)=nan;
    Wb(Wb==0)=nan;
    
    trace1_sim=Wb;
    trace1_real=Beta_hat;
    
    trace1_sim(isnan(trace1_sim))=0;
    trace1_real(isnan(trace1_real))=0;

    trace1_real(1)=0;
    trace1_sim(:,1)=0;
    
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
    eff_beta_hat1=eff_beta_hat((k-1)*n0+1:k*n0)';
    eff_beta_hat1(v<a|v>b)=0;
    Beta_hat0=[tril(ones(n0))*eff_beta_hat1*(v(3)-v(2))];
    fir=find(v>=a, 1, 'first'); las=find(v<=b, 1, 'last');
    Beta_hat1=(Beta_hat0-Beta_hat0(fir))./(v-v(fir))'-(Beta_hat0(las)-Beta_hat0(fir))/(v(las)-v(fir));%a cant be grid point
    Beta_hat2=sqrt(n)*Beta_hat1;
    Beta_hat=Beta_hat2.*(v>aa&v<b)';   %Q^{(1)}
    Beta_hat(isnan(Beta_hat))=0;
    
    vv=repmat(v,D,1);
    Wb0=G*hat_H(:,:,k)/sqrt(n);
    Wb1=(Wb0-repmat(Wb0(:,fir),1,n0))./repmat((v-v(fir)),D,1)-repmat((Wb0(:,las)-Wb0(:,fir))/(v(las)-v(fir)),1,n0);
    Wb=Wb1.*(vv>=aa&vv<=b);
    Wb(isnan(Wb))=0;
    cov_Wb=[0,cov_Wb1(:,k)'];    
    
    Beta_hat(Beta_hat==0)=nan;
    Wb(Wb==0)=nan;
        
    trace2_sim=Wb;
    trace2_real=Beta_hat;    
    
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
     
    M4_hat=nansum(-Beta_hat'.*(cov_Wb(2:n0+1)-cov_Wb(1:n0)));
    M4_tuta=nansum(-Wb.*repmat((cov_Wb(2:n0+1)-cov_Wb(1:n0)),D,1),2);
 
    pv_a1_t2(k)=mean(A1_tuta>A1_hat);
    pv_a2_t2(k)=mean(A2_tuta>A2_hat);
    pv_m1_t2(k)=mean(M1_tuta<M1_hat);
    pv_m2_t2(k)=mean(M2_tuta<M2_hat);
    pv_m3_t2(k)=mean(M3_tuta>M3_hat);
    pv_m4_t2(k)=mean(M2_tuta>M2_hat);
    
end
 


p_values_t1=[pv_a1;pv_a2;pv_m1;pv_m2;pv_m3;pv_m4];
p_values_t2=[pv_a1_t2;pv_a2_t2;pv_m1_t2;pv_m2_t2;pv_m3_t2;pv_m4_t2];

%%%%%%%%%%%%%%%% Cumulative incidence function  %%%%%%%%%%%%%%%%%%

if  (nstrt~=1 | Z_pd==99)  F=NaN;
else
    for k=1:size(Z,2)
        rr=ksrmv(v',eff_beta_hat((k-1)*length(v)+(1:length(v)))',v(3)-v(1),V);
        eff_beta_u(:,k)=rr.f;
    end
    eb=0;predeb=0;
    for k=1:size(Z,2)
        eb=repmat(eff_beta_u(:,k).*Z(:,k),1,n)+eb;
        predeb=repmat(eff_beta_u(:,k).*predZ(:,k),1,n)+predeb;
    end    
    
    n2=1000;
    v2=linspace(0.1,0.9,n2);
    d_v=diff(v2);    
    for j=1:n2
        for kk=1:nstrt
            S0kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n)))./nansum(repmat(strata==kk-1,1,n));
            dlambdakk=1./S0kk./sum(repmat(strata==kk-1,1,n));
            A0_tx(:,j)=nansum((repmat(X,1,n)>=repmat(X',n,1)).*repmat(delta',n,1).* repmat( K((V'-v2(j))/h)/h .*dlambdakk,n,1 ) ,2);
            %Var2kk(kk,:)=A0_tx'.*(strata==kk-1)' ;
        end
    end
    
    for j=1:n0
        eb0=0;
        for k=1:beta_n
            eb0=eff_beta_hat(j+n0*(k-1))*Z_pd(k)+eb0;
        end
        e_A=sum(A0_tx.*d_v(1)*exp(eb0),2);
        for kk=1:nstrt
            S0kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+sum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n)))./nansum(repmat(strata==kk-1,1,n));
            dlambdakk=1./S0kk./sum(repmat(strata==kk-1,1,n));      
            
            F(j)= sum( (X<=tt).*delta.*exp(-e_A).* K((V-v(j))/h)/h .*dlambdakk' )*exp(eb0) ;
        end
    end
    
end

%
%
%
% for k=1:size(Z,2)
%     rr=ksrmv(v',eff_beta_hat((k-1)*length(v)+(1:length(v)))',v(3)-v(1),V);
%     eff_beta_u(:,k)=rr.f;
% end
% eb=0;predeb=0;
% for k=1:size(Z,2)
%     eb=repmat(eff_beta_u(:,k).*Z(:,k),1,n)+eb;
%     predeb=repmat(eff_beta_u(:,k).*predZ(:,k),1,n)+predeb;
% end
%
% S0(1,:)=nanmean((repmat(X,1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+nanmean((repmat(X,1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n));
%
% dlambda=1/n./S0(1,:).*delta';
% dlambda(dlambda==0)=NaN;
% VV=repmat(v,n0,1);
%
% xx=linspace(0.1,max(X),n0);
%
% XX=repmat(xx',1,n0);
% lambda=XX;
% for iii=1:n0
%     for jjj=1:n0
%         lambda(iii,jjj)=nansum(K((xx(iii)-X)/h2)/h2.*K((v(jjj)-V)/h)/h.*dlambda');
%     end
% end
%
% lambda1=lambda.*(VV(1,2)-VV(1,1)).*(XX(2,1)-XX(1,1));
% for iii=1:n0
%     for jjj=1:n0
%         Lambda(iii,jjj)=sum(sum(lambda1(1:iii,1:jjj)));
%     end
% end
%
%
% xx=linspace(0.1,max(X),n0);
% h2=0.4;
%
% for kk=1:nstrt
%     S0kk=(nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(eb).*repmat(wei,1,n))+nansum((repmat(X.*(strata==kk-1),1,n)>=repmat(X',n,1)).*exp(predeb).*repmat(1-wei,1,n)))./sum(repmat(strata==kk-1,1,n));
%     S0(kk,:)=S0kk.*(strata==kk-1)';
%
%     dlambda=1/n./S0(kk,:).*delta';
%     dlambda(dlambda==0)=NaN;
%     VV=repmat(v,n0,1);
%     XX=repmat(xx',1,n0);
%     lambda=XX;
%     for iii=1:n0
%         for jjj=1:n0
%             lambda(iii,jjj)=nansum(K((xx(iii)-X)/h2)/h2.*K((v(jjj)-V)/h)/h.*dlambda');
%         end
%     end
%
%     lambda1=lambda.*(VV(1,2)-VV(1,1)).*(XX(2,1)-XX(1,1));
%     for iii=1:n0
%         for jjj=1:n0
%             Lambda(iii,jjj,kk)=sum(sum(lambda1(1:iii,1:jjj)));
%         end
%     end
%
% end
%
%
