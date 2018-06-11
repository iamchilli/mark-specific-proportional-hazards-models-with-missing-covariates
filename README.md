# mark-specific-proportional-hazards-models-with-missing-covariates
Hypothesis tests for stratified mark‐specific proportional hazards models with missing covariates, with application to HIV vaccine efficacy trials

## Instructions for using the MatLab code for simulations 


The following are the instructions for using the matlab program code. The tables for the simulation results in the manuscript were produced without using a fixed set of random seeds. Because of this, it will not be able to produce the exact simulation results. However, the results are reproducible (i.e. consistent) to the extend of the differences caused by random seeds. The required computation time for each setting for sample size 500 is about 2 hours running on the High Performance Computing Cluster at UNC Charlotte.

The scripts used to produce Tables 1-8 are named table1.m to table8.m, respectively. The core function used in all these scripts is main.m function. The MatLab files aipw_strt.m, aipw.m, cuminc.m, ipw_strt.m, ksrmv.m, simulate_2noT.m, simulate_3noT.m, simulate.m and weight.m are the auxilliary functions.

The input of the main.m function is:

    % nsim   Number of simulations

    % n       Sample size

    % n0                 % Number of grid points

    % gamma                 % gamma cannot be zero, if strata>1, it could be gamma=[-0.6, -0.3];

    % alpha              % Model parameter

    % beta1              % Model parameter

    % beta2              % Model parameter

    % beta3              % Model parameter, should be [] for table 1-6

    % n_strat            % Number of stata, or n_strat=1, 2 or more

    % strata=binornd(1,0,n,1);% strata levels: 0, 1, 2 ..... (Note: It must contain zero)

    % num_cov=3          % Number of parameters

    % aux_corr   	     %Correlation with aux variables 0, 0.7 or 0.9 

    % model              % Model for generate phase 2 indicator


The output  of the main.m function is:

    %plot_aipw_beta		%estimation

    %power_a1 		%power of test 1

    %power_t2_a1		%power of test 2



To use the code to produce the simulation results, please start with the tableX.m files as follows:

    [plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, n, n0, gamma, alpha , beta1 , beta2 , beta3 , n_strat ,strata, num_cov,aux_corr, model)


We have provided the intermediate results that were used for Tables 1-8. They are named as tableX.mat (table1a.mat, table1b.mat, etc.). The following code could get the results for Tables 1 of the manuscript.

    load table1a.mat ;
    fprintf('%12.3f', out )
    fprintf('\n'   )
    load table1b.mat ;
    fprintf('%12.3f', out )
    fprintf('\n'   )
    load table1c.mat ;
    fprintf('%12.3f', out )
    fprintf('\n'   )
