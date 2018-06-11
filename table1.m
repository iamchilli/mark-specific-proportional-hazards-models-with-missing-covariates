
%Each block can be runned parallel. 

clear all
nsim=500
%%%%%%%%%%%%%%% table1 %%%%%%%%%%%%%%%%%%%%
%table 1, rho=0, M1, n=500
[plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, 500, 20, 0.6, 0, 0, -0.3,  [], 1, binornd(1,0,500,1), 2, 0, 1);
out1=[ power_a1(1), power_a2(1), power_m1(1), power_m2(1)];
fprintf('%12.3f',  out1)
 

%table 1, rho=0, M1, n=800
[plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, 800, 20, 0.6, 0, 0, -0.3, [], 1, binornd(1,0,800,1), 2, 0, 1);
out2=[   power_a1(1), power_a2(1), power_m1(1), power_m2(1)];
fprintf('%12.3f',  out2)
 

%table 1, rho=0, M2, n=500
[plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, 500, 20, 0.6, -0.7, 0.7, -0.3, [], 1, binornd(1,0,500,1), 2, 0, 1);
out3=[   power_a1(1), power_a2(1), power_m1(1), power_m2(1)];
fprintf('%12.3f',  out3)
 

%table 1, rho=0, M2, n=800
[plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, 800, 20, 0.6, -0.7, 0.7, -0.3, [], 1, binornd(1,0,800,1), 2, 0, 1);
out4=[   power_a1(1), power_a2(1), power_m1(1), power_m2(1)];
fprintf('%12.3f',  out4)
 

%table 1, rho=0, M3, n=500
[plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, 500, 20, 0.6, -1.5, 1.5, -0.3, [],  1, binornd(1,0,500,1), 2, 0, 1);
out5=[   power_a1(1), power_a2(1), power_m1(1), power_m2(1)];
fprintf('%12.3f',  out5)
 

%table 1, rho=0, M3, n=800
[plot_aipw_beta,power_a1,power_a2,power_m1,power_m2,power_t2_a1,power_t2_a2,power_t2_m1,power_t2_m2]=main(nsim, 800, 20, 0.6, -1.5, 1.5, -0.5, [], 1, binornd(1,0,800,1), 2, 0, 1);
out6=[   power_a1(1), power_a2(1), power_m1(1), power_m2(1)];
fprintf('%12.3f',  out6)

out=[out1,out2,out3,out4,out5,out6]
fprintf('%12.3f',  out)
save table1a.mat 

%%%% The rest of the blocks of table 1 can be got by changing aux_corr=0.7 and aux_corr=0.93 (please check readme file)%%%%
