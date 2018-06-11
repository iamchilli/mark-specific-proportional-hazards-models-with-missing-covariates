function [plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim,n, n0, gamma, alpha , beta1 , beta2 , beta3 , n_strat ,strata, num_cov,aux_corr,model)

%progress bar 
% bar = waitbar(0,'1','Name','Number of Simulations done',...
%             'CreateCancelBtn',...
%             'setappdata(gcbf,''canceling'',1)');
% setappdata(bar,'canceling',0)

% nsim=5;               % Number of simulations
% n=500;                  % Sample size
% n0=33;                  % Number of grid points
% 
% gamma=[0.3];           % gamma cannot be zero, if strata>1, it could be gamma=[-0.6, -0.3];
% alpha=0;              % Model parameter
% beta1=0;              % Model parameter
% beta2=-0.3;             % Model parameter
% beta3=0.5;              % Model parameter
% n_strat=1;              % or n_strat=2 or more
% strata=binornd(1,0,n,1);% strata levels: 0, 1, 2 ..... (Note: It must contain zero)
% num_cov=3;



v=linspace(0.1,0.9,n0); % Interval of grid points of v, for example ,[0.1,0.9];

for i=1:nsim
    %progress bar 
    %  if getappdata(bar,'canceling')
    %     break
    % end
    % % Report current estimate in the waitbar's message field
    % waitbar(i / nsim,bar,sprintf('%12.0f',i))
 
    
    rng(2*i+1)
    %[Z,Zt,Aux,X,V,delta]=simulate(n,alpha,gamma,beta1,beta2, beta3, strata);    % 3 variables with one time depdent variable
   if (num_cov==3)
        [Z,Zt,Aux,X,V,delta]=simulate_3noT(n,alpha,gamma,beta1,beta2, beta3,strata); % 3 variables with no time depdent variable
   end
   if (num_cov==2)
        [Z,Aux,X,V,delta ]=simulate_2noT(n,alpha,gamma,beta1,beta2,aux_corr);
        Zt=[];
    end    

    delta=delta+0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Part II: Method specication       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%% bandwidth and Kernel %%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_op = 4*nanstd(V) * length(V(~isnan(V),:))^(-1/3);   % output bandwidth by formula
    h=0.15;                                           % bandwidth to be used
    K=@(x)0.75*(1-x.^2).*(abs(x)<=1);
    
    %%%%%%%%%%%%%%%%%% missing covariate missing pattern %%%%%%%%%%
    delta=delta+0;
    if (model==1) sel=binornd(1,0.7,n,1);  end % MAR with p=0.7
    if (model==2) sel=binornd(1, delta*0.8+(1-delta)*0.6,n,1); end 
    if (model==3) sel=binornd(1,exp(1.2*Z(:,1)-0.4*delta-0.3)./(1+exp(1.2*Z(:,1)-0.4*delta-0.3)),n,1); end
    %sel=binornd(1,exp(Z(:,2)+delta-0.5)./(1+exp(Z(:,2)+delta-0.5)),n,1);
    
    
    %%%%%%%%%%%%%%%%%% Interval limits used in the tesing procedures %%%%%%%%%%
    ll1=0.11;
    ul=0.89;
    ll2=0.21;
    
    %%%%%%%%%%%%%%%%%% Used for calculation cumulative incidence %%%%%%%%%%
    if (num_cov==3) Z_pd=[0.5,0.5,0.5]; end % used for calculation cumulative incidence
    if (num_cov==2) Z_pd=[0.5,0.5]; end % used for calculation cumulative incidence
    tt=1;                % used for calculation cumulative incidence
    
    beta_Z=size(Z,2);
    beta_Zt=size(Zt,2)/n;
    
    %%%%%%%%%%%%%%%%%% Weight and conditioanl expection of the AIPW estimation %%%%%%%%%%
    method='parametric';
    Z_c=[Z(:,2) delta Aux  ];
    [wei]=weight([Z(:,2) delta Aux],sel,method);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Part III: IPW and AIPW estimation  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IPW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %beta_ini=repmat(0,beta_Z+beta_Zt,size(v,2));%initial value
    %[beta_hat_ipw,sig_beta_ipw]=ipw_strt(beta_ini,X,Z,Zt,V,v,delta,sel,beta_Z,beta_Zt,wei,h,n0,n,K,strata);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AIPW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta_ini=repmat(0,beta_Z+beta_Zt,size(v,2));%initial value
    [beta_hat_aipw,sig_beta_aipw,p_values_t1,p_values_t2,trace1_sim,trace1_real,trace2_sim,trace2_real,F]=aipw_strt(beta_ini,X,Z,Zt,V,v,delta,sel,beta_Z,beta_Zt,wei,h,n0,n,K,Z_c,method,ll1,ul,ll2,strata, Z_pd, tt);
    
    aipw_beta_hat(i,:)=beta_hat_aipw;
    sig_aipw_beta(i,:)=sig_beta_aipw;
    
    p_value_a1(i,:)=p_values_t1(1,:);
    p_value_a2(i,:)=p_values_t1(2,:);
    p_value_m1(i,:)=p_values_t1(3,:);
    p_value_m2(i,:)=p_values_t1(4,:);
    
    
    p_value_t2_a1(i,:)=p_values_t2(1,:);
    p_value_t2_a2(i,:)=p_values_t2(2,:);
    p_value_t2_m1(i,:)=p_values_t2(3,:);
    p_value_t2_m2(i,:)=p_values_t2(4,:);
  
end

if (num_cov==3)

beta10=beta1*v.*v+alpha;
beta20=beta2*v.*v;
beta_true=[beta10,beta20,beta3*(v>0)];
end

if (num_cov==2)

beta10=beta1*v+alpha;
beta20=beta2*v;
beta_true=[beta10,beta20 ];
end


EST_aipw_beta=mean(aipw_beta_hat);
SE_aipw_beta=mean(sig_aipw_beta);
SD_aipw_beta=std(aipw_beta_hat);
CP_aipw_beta=nanmean(aipw_beta_hat>=repmat(beta_true,nsim,1)-1.96*sig_aipw_beta & aipw_beta_hat<=repmat(beta_true,nsim,1)+1.96*sig_aipw_beta);
plot_aipw_beta=[EST_aipw_beta-beta_true;SD_aipw_beta;SE_aipw_beta;CP_aipw_beta];



power_a1=1-mean(p_value_a1>=0.05);
power_a2=1-mean(p_value_a2>=0.05);
power_m1=1-mean(p_value_m1>=0.05);
power_m2=1-mean(p_value_m2>=0.05);

power_t2_a1=1-mean(p_value_t2_a1>=0.05);
power_t2_a2=1-mean(p_value_t2_a2>=0.05);
power_t2_m1=1-mean(p_value_t2_m1>=0.05);
power_t2_m2=1-mean(p_value_t2_m2>=0.05);




% delete(bar)
