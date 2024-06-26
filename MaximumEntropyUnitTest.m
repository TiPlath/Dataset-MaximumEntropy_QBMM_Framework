clc
clear variables
close all

%% Maximum Entropy reconstruction unit test
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Goal of this unit test is to plot the maximum entropy reconstruction with
% the same function it uses to reconstruct. We first assemble a function,
% compute the moments and then try to get the same moments + reconstruction
% back. That will show us if the maximum entropy reconstruction works
% correctly

%% Compute MER for 3 lambda assuming this is a volume density distribution (v_L)
expFunction2L = @(lambda,x) exp(-(lambda(1)+lambda(2)*x+lambda(3)*x.^2));
expFunction2LScaled =  @(lambda,x,M0) exp(-(lambda(1)+lambda(2)*x/max(x) + lambda(3)*(x/max(x)).^2))*M0/max(x);
lambda = [1,1,1];
% number of moments
N = length(lambda);
x = linspace(0.01,10,300);
y = expFunction2L(lambda,x);
% plot function
figure(1)
plot(x,y,"DisplayName","analytic2Lambda",LineWidth=2)
hold on

M = ComputeMoments(x,y,N);

[PSD,PSDErrorComp,k,E,lambdaME] = getPSD(x,M,N);

plot(x,PSDErrorComp(k,:),"--","DisplayName","ME2Lambda",LineWidth=2,MarkerSize=6)

%% reassemble function with LambdaME

yME = expFunction2LScaled(lambdaME,x,M(1));

MME = ComputeMoments(x,yME,N);

plot(x(1:10:end),yME(1:10:end),"o","DisplayName","ME2LambdaReassembled",LineWidth=2,MarkerSize=6)

legend()

%% Doing the same with a number density distribution n_L (DOES NOT WORK WELL)
% volume = x;
% nL_expFunction2L = @(lambda,x) (exp(-(lambda(1)+lambda(2)*x))) ./ volume;
% nL_expFunction2LScaled =  @(lambda,x,M0) (exp(-(lambda(1)+lambda(2)*x/max(x)))*M0/max(x));
% 
% nL_y = nL_expFunction2L(lambda,x);
% %nL_y(1) = 2*nL_y(2);
% % plot function
% figure(2)
% plot(x,nL_y,"DisplayName","analytic2Lambda n_L",LineWidth=2)
% hold on
% 
% nL_M = ComputeMoments(x,nL_y,N);
% 
% [nL_PSD,nL_PSDErrorComp,nL_k,nL_E,nL_lambdaME] = getPSD(x,nL_M,N);
% 
% plot(x,nL_PSDErrorComp(k,:),"--","DisplayName","ME2Lambda",LineWidth=2,MarkerSize=6)
% 
% plot(x,PSDErrorComp(k,:) ./ volume, "-x","DisplayName","ME2Lambda v_L to n_L",LineWidth=1,MarkerSize=10)
% 
% plot(x(1:end), yME(1:end) ./ volume,"-o","DisplayName","ME2LambdaReassembled v_L to n_L",LineWidth=1,MarkerSize=4)
% 
% %% reassemble function with LambdaME
% 
% nL_yME = nL_expFunction2LScaled(nL_lambdaME,x,nL_M(1));
% 
% nL_MME = ComputeMoments(x,nL_yME,N);
% 
% plot(x(1:end),nL_yME(1:end),"o","DisplayName","ME2LambdaReassembled",LineWidth=2,MarkerSize=6)
% % plot(x,y,"x","DisplayName","analytic2Lambda v_L",LineWidth=2)
% 
% legend()
% 
% %% add reconstructed figure(2) distributions (n_L) as v_L distribution to figure(1) (THAT DOES NOT WORK! DO NOT RECONSTRUCT AND SWITCH DISTRIBUTION AFTERWARDS!)
% 
% figure(1)
% hold on
% 
% plot(x,nL_PSD(k,:).*volume,"-x","DisplayName","ME2Lambda n_L to v_L",LineWidth=1,MarkerSize=4)
% plot(x(1:end),nL_yME(1:end).*volume,"-o","DisplayName","ME2LambdaReassembled n_L to v_L",LineWidth=1,MarkerSize=4)
