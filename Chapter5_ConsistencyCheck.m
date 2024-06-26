%% Application to twin-screw wet granulation data: Consistency check (2) (section 5.3)

% We perform a consistency check our created model form the dataset analysis 
% (1) using the twin-screw wet granulation dataset (Plath et al. 2021) with
% the lactose, polyvinylpyrrolidon and microcrystalline cellulose mixture.
% 1. For one of the three centerpoint short-screw experiments (N28), we 
%    predict the mean residence time MRT as well as the functions found in the
%    dataset analysis (1) to predict the Sauter mean diameter d_32, minimum 
%    aggregation rate a_s and non-dimensional breakage rate beta_t'. 
% 2. We use the same kernel values to predict the centerpoint long-screw 
%    experiment (N22) by simulating for another mean residence time, 
%    confirming steady state. The reconstructed particle size distribution
%    is compared to the experimental distributions to check consistency of 
%    our DQMOM-MER framework.

clc
close all
clear variables
% disable legend warnings
warningState = warning('off', 'MATLAB:legend:IgnoringExtraEntries');
% disable nearly singular matrix warning
warning('off', 'MATLAB:nearlySingularMatrix');
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define important variables (pre-processing)
% The total volume [m^3] of a typical sample from twin-screw wet granulation
V_P = 2.05e-6;
% Total number of particles (zeroth moment)
N_P = 8.78888e+07;
% Density
rho_P = 1540;
% number of times to reconstruct
ntimes = 2;
% resolution of the reconstructed PSD by MER
resolution = 300;
% time step
dtPrime=0.001;
% Number of dirac-delta distributed classes (weights and nodes, 1,...,7)
N_delta=5;
% Maximum number of moments
mMax = 2*N_delta;
% define process parameter pair
SFL = 0.08;
LS = 0.225;
experimentFit = "MCCLactosePVP_CenterPointShortScrew_N28_Validation";
% Number of fragments per broken particle length-based
N_fL =8;
%% read custom Distribution
% call function to read distribution (v_L CDF)
CDFv_L = ReadData('Data/InitialLactoseMCCPVP_CDF-v_L.csv',1e6,0,0.995);
CDFv_L_Fit = ReadData("Data/" + experimentFit + ".csv",1e6,0,1);
CDFv_L_Predict = ReadData('Data/MCCLactosePVP_CenterPointLongScrew_N22_Validation.csv',1e6,0,1);
% Normalize and scale by V_P to represent volume of a typical wet
% granulation sample (assume volume should stay constant)
CDFv_L(:,2) = CDFv_L(:,2)/max(CDFv_L(:,2));
CDFv_L(:,2) = CDFv_L(:,2)*V_P;
CDFv_L_Fit(:,2) = CDFv_L_Fit(:,2)/max(CDFv_L_Fit(:,2));
CDFv_L_Fit(:,2) = CDFv_L_Fit(:,2)*V_P;
CDFv_L_Predict(:,2) = CDFv_L_Predict(:,2)/max(CDFv_L_Predict(:,2));
CDFv_L_Predict(:,2) = CDFv_L_Predict(:,2)*V_P;
%% Density distribution conversions
% Convert CDF v_L to PDF v_L for the MER
v_L = convertCDFtoPDF(CDFv_L);
v_L_Fit = convertCDFtoPDF(CDFv_L_Fit);
v_L_Predict = convertCDFtoPDF(CDFv_L_Predict);
% convert CDFv_L to PDF n_V for the DQMOM
% Convert CDF v_L to volume-based VDF (v_V)
kappa = pi/6;
CDFv_V(:,2) = CDFv_L(:,2);
CDFv_V(:,1) = kappa .* CDFv_L(:,1).^3;
CDFv_V_Fit(:,2) = CDFv_L_Fit(:,2);
CDFv_V_Fit(:,1) = kappa .* CDFv_L_Fit(:,1).^3;
CDFv_V_Predict(:,2) = CDFv_L_Predict(:,2);
CDFv_V_Predict(:,1) = kappa .* CDFv_L_Predict(:,1).^3;
% Convert CDF v_V to PDF v_V
v_V = convertCDFtoPDF(CDFv_V);
v_V_Fit = convertCDFtoPDF(CDFv_V_Fit);
v_V_Predict = convertCDFtoPDF(CDFv_V_Predict);
% Convert volume-based VDF (v_V) to volume-based NDF (n_V)
n_V(:,1) = v_V(:,1);
n_V(:,2) = v_V(:,2)./v_V(:,1);
n_V_Fit(:,1) = v_V_Fit(:,1);
n_V_Fit(:,2) = v_V_Fit(:,2)./v_V_Fit(:,1);
n_V_Predict(:,1) = v_V_Predict(:,1);
n_V_Predict(:,2) = v_V_Predict(:,2)./v_V_Predict(:,1);
% Conversion to get n_L
% convert PDF n_V to CDF n_V
CDFn_V = convertPDFtoCDF(n_V);
CDFn_V_Fit = convertPDFtoCDF(n_V_Fit);
CDFn_V_Predict = convertPDFtoCDF(n_V_Predict);
% convert CDFn_V to CDFn_L
CDFn_L = CDFn_V;
CDFn_L(:,1) = (CDFn_V(:,1)/kappa).^(1/3);
CDFn_L_Fit = CDFn_V_Fit;
CDFn_L_Fit(:,1) = (CDFn_V_Fit(:,1)/kappa).^(1/3);
CDFn_L_Predict = CDFn_V_Predict;
CDFn_L_Predict(:,1) = (CDFn_V_Predict(:,1)/kappa).^(1/3);
% convert CDFn_L to PDF n_L
n_L = convertCDFtoPDF(CDFn_L);
n_L_Fit = convertCDFtoPDF(CDFn_L_Fit);
n_L_Predict = convertCDFtoPDF(CDFn_L_Predict);
%% Define Scale variables for non-dimensionalisation at t=0 (pre-processing)
% Compute total number of particles to define particle unit scale P
M_0 = sum(n_V(1:end-1,2).*diff(n_V(:,1)));
% Compute total volume to define length scale L
M_1V = sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1)));
% Particle scale [p]
Omega_P = M_0;
M0Prime = 1;
% Volume scale [m^3]
Omega_V = (M_1V * Omega_P^(-1));
M1Prime = 1;

%% Predictions for MRT, d_32, a_s and beta_tPrime found from dataset analysis
% mean residence time prediction from SFL and LS (short-screw)
MRTFunc = @(params, xy) params(1) * xy(:,1)+ params(2) * xy(:,2) + params(3);
abc1 = [123.46,45.16,3.89];
parameterPair = [SFL, LS];
MRT = MRTFunc(abc1,parameterPair);

% Sauter mean diameter prediction from SFL and LS (short-screw)
d32Func = @(params,xy) (params(1).*xy(:,2).^2) .* (xy(:,1)./(1000*rho_P)).^(1/3);
abc2 = 6.06368;
d32_fitFunction = d32Func(abc2,parameterPair);

% minimum aggregation rate prediction from MRT, d32
t_opt = 7.5;
a_sFunc = @(params, xy) params(1) * ((xy(:,1)-t_opt)/t_opt).^(params(2)) .* (xy(:,2)*Omega_V^(-1/3)).^(params(3));
% a_c, p_t, p_m
abc3 = [0.3184, -0.55, -0.181];
a_s = a_sFunc(abc3, [MRT,d32_fitFunction]);  
a_t = a_s*(Omega_P*Omega_V)^(-1);

% non-dimensional breakage rate prediction from d32
beta_tPrimeFunc = @(params, x) params(1).*x.^(-params(2));
abc4 = [2.625, 0.914];
d32Prime_fitFunction = d32_fitFunction*Omega_V^(-1/3);
beta_tPrime = beta_tPrimeFunc(abc4,d32Prime_fitFunction);
%% set rates for the kernels using the predictions and compute timescale
% Aggregation rate
Omega_a = a_t*Omega_V;  
aPrime(1:N_delta,1) = 1;
% Time scale, where aggregation is assumed to be dominant [s]
Omega_t = Omega_P^(-1)*Omega_a^(-1);
% Growth rate length-based
gPrime_L(1:N_delta,1) = 0.00152;
% Breakage rate length-based
betaPrime_L(1:N_delta,1) = beta_tPrime;
% symmetric fragentation for a length-based PSD (here: n_L)
b_alpha_L = @(L_alpha,k,N_fL) N_fL^((3-k)/3)*L_alpha.^k;

%% determine maximum simulation time and initialise parameters
% max time is similar to 2x MRT
tmaxPrime=round(2*MRT*Omega_t^(-1),2);
% Number of time steps
ntPrime=round(tmaxPrime/dtPrime)+1;
% Initialize different iteration vectors (useful for debugging and plotting)
% A matrix storing the volume-based momenta of every time step
iterM_L=zeros(ntPrime,2*N_delta);
% store weights and nodes of every iteration
iterL = zeros(ntPrime,N_delta);
iterw_L = zeros(ntPrime,N_delta);
% Initialize x-vector to store domain data of every iteration
iterx = zeros(ntimes+1,resolution);
% store PSD data of every iteration
iterPSD = zeros(ntimes+1,resolution);
% store lambdas of every iteration
iterLambda = zeros(ntimes+1,2*N_delta);
% store dxMax of every iteration
iterdxMax = zeros(ntPrime,1);

%% Define scaled variables
% define dimensionless particle Volume V [m^3] and distribution function n_V [p/m^3]
% to ensure we get dimensionless moments M_V
n_VPrime(:,1) = n_V(:,1)./(Omega_V);
n_VPrime(:,2) = n_V(:,2).*Omega_V.*Omega_P^(-1);
n_VPrime_Fit(:,1) = n_V_Fit(:,1)./(Omega_V);
n_VPrime_Fit(:,2) = n_V_Fit(:,2).*Omega_V.*Omega_P^(-1);
n_VPrime_Predict(:,1) = n_V_Predict(:,1)./(Omega_V);
n_VPrime_Predict(:,2) = n_V_Predict(:,2).*Omega_V.*Omega_P^(-1);
% define dimensionless particle Length L [m] and distribution function v_L [p m^3/m]
% to ensure we get dimensionless moments Mv_L
v_LPrime(:,1) = v_L(:,1)./(Omega_V)^(1/3);
v_LPrime(:,2) = v_L(:,2).*Omega_V^(-2/3).*Omega_P^(-1);
v_LPrime_Fit(:,1) = v_L_Fit(:,1)./(Omega_V)^(1/3);
v_LPrime_Fit(:,2) = v_L_Fit(:,2).*Omega_V^(-2/3).*Omega_P^(-1);
v_LPrime_Predict(:,1) = v_L_Predict(:,1)./(Omega_V)^(1/3);
v_LPrime_Predict(:,2) = v_L_Predict(:,2).*Omega_V^(-2/3).*Omega_P^(-1);
% define dimensionless n_L
n_LPrime(:,1) = n_L(:,1)./(Omega_V^(1/3));
n_LPrime(:,2) = n_L(:,2).*Omega_V^(1/3).*Omega_P^(-1);
n_LPrime_Fit(:,1) = n_L_Fit(:,1)./(Omega_V^(1/3));
n_LPrime_Fit(:,2) = n_L_Fit(:,2).*Omega_V^(1/3).*Omega_P^(-1);
n_LPrime_Predict(:,1) = n_L_Predict(:,1)./(Omega_V^(1/3));
n_LPrime_Predict(:,2) = n_L_Predict(:,2).*Omega_V^(1/3).*Omega_P^(-1);
% dimensionless time
tPrime=(0:ntPrime-1)*dtPrime;

%% get initial weights and nodes for a volume-based distribution (v_L)
% Compute length-based moments of distribution for the initial reconstruction
Mv_LPrime = ComputeMoments(v_LPrime(:,1), v_LPrime(:,2), mMax);
Mv_LPrime_Fit = ComputeMoments(v_LPrime_Fit(:,1), v_LPrime_Fit(:,2), mMax);
Mv_LPrime_Predict = ComputeMoments(v_LPrime_Predict(:,1), v_LPrime_Predict(:,2), mMax);
% we compute weights and nodes using the Wheeler algorithm
[Lv_LPrime,wv_LPrime] = Wheeler(Mv_LPrime(1:2*N_delta),N_delta);
[Lv_LPrime_Fit,wv_LPrime_Fit] = Wheeler(Mv_LPrime_Fit(1:2*N_delta),N_delta);
[Lv_LPrime_Predict,wv_LPrime_Predict] = Wheeler(Mv_LPrime_Predict(1:2*N_delta),N_delta);
%% reconstruct the initial density distribution by maximum entropy approach using v_L
[d,initE,M_,initx,initPSD,c,lambda] = plotInitLeastErrorPSD(v_LPrime(:,1),v_LPrime(:,2),Lv_LPrime,Mv_LPrime,mMax);
%% plot final distributions
%close figures to override ploting options
clf
figure(1)
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log','FontSize',16);
set(gcf, 'Color', 'w');
xlabel('$L$','FontSize',18,'Interpreter','Latex')
ylabel('$v_L$','FontSize',18,'Interpreter','Latex')
hold on
plot(v_LPrime_Fit(:,1),v_LPrime_Fit(:,2), 'LineWidth', 2);
plot(v_LPrime_Predict(:,1),v_LPrime_Predict(:,2), 'LineWidth', 2);
% fix xMin to the minimum particle size of the initial reconstruction
xMin = 0;
% store initial densities, reconstruction domain to the iteration vectors
iterx(1,1:end) = initx;
iterPSD(1,1:end) = initPSD;
iterLambda(1,1:length(lambda)) = lambda;
% Legend for plotting
str = {'$v_{\mathrm{L}}^{\prime}$ short-screw', '$v_{\mathrm{L}}^{\prime}$ long-screw'};

%% evolve abscissae and weights in time according to DQMOM using n_L
% Compute length-based moments of distribution n_L
Mn_LPrime = ComputeMoments(n_LPrime(:,1), n_LPrime(:,2), mMax);
Mn_LPrime_Fit = ComputeMoments(n_LPrime_Fit(:,1), n_LPrime_Fit(:,2), mMax);
Mn_LPrime_Predict = ComputeMoments(n_LPrime_Predict(:,1), n_LPrime_Predict(:,2), mMax);
% We compute weights and nodes from n_L using the Wheeler algorithm
[Ln_LPrime,wn_LPrime] = Wheeler(Mn_LPrime(1:2*N_delta),N_delta);
% Store inital moments, nodes and weights of the Volume based distribution for
% further processing (n_L)
iterM_L(1,:) = getMomenta(Ln_LPrime,wn_LPrime);
% Store inital nodes and weights in iteration vectors
iterL(1,1:N_delta) = Ln_LPrime;
iterw_L(1,1:N_delta) = wn_LPrime;
iterdxMax(1) = abs(max(Ln_LPrime)*(c-1));
iterd32(1) = Mn_LPrime(4)/Mn_LPrime(3);
% define counter
o = 0;
% Time loop, forward Euler
for i=2:ntPrime
    % Compute one time step of DQMOM (n_L)
    [Ln_LPrime,wn_LPrime] = DirectQuadratureMethodOfMomentsNonDimensionalizedLengthBased(iterL(i-1,1:N_delta)',iterw_L(i-1,1:N_delta)',dtPrime,gPrime_L,aPrime,betaPrime_L,b_alpha_L,N_fL);
    % compute momenta from weights and nodes (n_L)
    iterM_L(i,:) = getMomenta(Ln_LPrime,wn_LPrime);
    % Assign weights and nodes to iteration vectors
    iterL(i,1:N_delta) = Ln_LPrime;
    iterw_L(i,1:N_delta) = wn_LPrime;
    % compute distance to extend over the actual function and save it
    dxMax = abs(max(Ln_LPrime)*(c-1));
    iterdxMax(i) = dxMax;
    if abs(tPrime(i) - round(MRT/Omega_t,2)) <= 1e-4 || abs(tPrime(i) - round(2*MRT/Omega_t,2)) <= 1e-4  
        o = o + 1;
        Lv_LPrime = (iterL(i,1:N_delta));
        wv_LPrime = iterw_L(i,1:N_delta).*iterL(i,1:N_delta).^3;
        % Reconstruct and plot the PSD from the distributions moments
        [PSD,x,k,lambda] = plotMaximumEntropyReconstruction(Lv_LPrime,wv_LPrime,mMax,xMin,dxMax,v_LPrime(:,2));
        xlabel('$L^{\prime}$','Interpreter','Latex')
        ylabel('$v_{\mathrm{L}}^{\prime}$','Interpreter','Latex')        
        if o == 1
            str = [str , sprintf('$\\bar{v}_{\\mathrm{L}}^{\\prime} (t^{\\prime} = \\mathit{MRT}_1)$' , tPrime(i))];
            iterPSD(2,1:length(PSD)) = PSD;
            iterx(2,1:length(x)) = x;
            iterLambda(2,1:length(lambda)) = lambda;
        end
        if o == 2
            str = [str , sprintf('$\\bar{v}_{\\mathrm{L}}^{\\prime} (t^{\\prime} = \\mathit{MRT}_2)$' , tPrime(i))];
            iterPSD(3,1:length(PSD)) = PSD;
            iterx(3,1:length(x)) = x;
            iterLambda(3,1:length(lambda)) = lambda;
        end
        legend(str, 'location', 'best','Interpreter','Latex')
        % Store number of moments used in its iteration vector
        iterk(i) = k;
    end    
end

%% Postprocessing 
%% plot d32 match
stepsize = 10;
d32Prime_Fit = (Mn_LPrime_Fit(4))/(Mn_LPrime_Fit(3));
d32Prime_predict = (Mn_LPrime_Predict(4))/(Mn_LPrime_Predict(3));
iterd32 = zeros(1,ntPrime);
iterd32 = (iterM_L(:,4))./(iterM_L(:,3));
d32_short =(iterM_L(floor(ntPrime/2),4))./(iterM_L(floor(ntPrime/2),3));
d32_long =(iterM_L(end,4))./(iterM_L(end,3));
figure(13)
axes13 = axes('Parent', figure(13));
plot(tPrime(1:stepsize:ntPrime),iterd32(1:stepsize:ntPrime), "x", 'MarkerSize', 14, 'LineWidth',2)
hold on
xlabel('$t^{\prime}$','Interpreter','Latex')
ylabel('$d_{32}^{\prime}$','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes13,'FontSize', 16)
line([axes13.XLim(1), axes13.XLim(2)], [d32Prime_Fit, d32Prime_Fit], 'Color', 'red', 'LineStyle', '-', 'LineWidth', 2);
legend('$d_{32}^{\prime}$', '$d_{32,\mathrm{Fit}}^{\prime}$', 'location','best', 'FontSize', 14,'Interpreter','Latex')

d32Error = (iterd32(end)-d32Prime_Fit)/d32Prime_Fit * 100;
fprintf("Sauter mean diameter is fit with an error of: %.2f%%\n\n", d32Error)
d32Prime = iterd32(end);
%% compute Jenson-Shannon divergence for simulated short- and long-screw distributions
%interpolate to Fit data to get same amount of points
v_LPrime_Fit_interp = interp1(v_LPrime_Fit(:,1),v_LPrime_Fit(:,2),iterx(2,:), 'linear', 0);
R_Fit = jensenShannonDivergence(iterPSD(2,:),v_LPrime_Fit_interp);
% R_Fit =  norm(iterPSD(2,:) - v_LPrime_Fit_interp)/iterM_L(1,1);
fprintf('L2-Norm for the fit reconstruction: R_Fit = %g\n', R_Fit)
%interpolate to Predict data to get same amount of points
v_LPrime_Predict_interp = interp1(v_LPrime_Predict(:,1),v_LPrime_Predict(:,2),iterx(3,:), 'nearest',0);
R_Predict = jensenShannonDivergence(iterPSD(2,:),v_LPrime_Predict_interp);
% R_Predict =  norm(iterPSD(3,:) - v_LPrime_Predict_interp)/iterM_L(1,1);
fprintf('L2-Norm for the predicted reconstruction: R_Predict = %g\n',R_Predict);

%% scale the moments to the same units
for k = 1:mMax-1
    iterM_L(:,k+1) = iterM_L(:,k+1).^(1/(k));
    Mn_LPrime_Fit(k+1) = Mn_LPrime_Fit(k+1).^(1/(k));
    Mn_LPrime_Predict(k+1) = Mn_LPrime_Predict(k+1).^(1/(k));
end
%% Plot moment error of short-screw distribution
for i = 0:mMax-1
    E_L_fit(i+1) = abs((iterM_L((ntPrime),i+1)/Mn_LPrime_Fit(i+1)) - 1);
end
fprintf('relative mean moment error for the d32 fit is: E = %g\n', sum(E_L_fit)/length(E_L_fit))

%% Plot moment error of long-screw distribution
for i = 0:mMax-1
    E_L_Predict(i+1) = abs((iterM_L((ntPrime),i+1)/Mn_LPrime_Predict(i+1)) - 1);
end
fprintf('relative mean moment error for the d32 prediction is: E = %g\n\n', sum(E_L_Predict)/length(E_L_Predict))

%% plot moments against moment order semilogarithmic
momentOrder = 0:1:(size(iterM_L,2) - 1);
figure(16)
axes16 = axes('Parent', figure(16));
semilogy(momentOrder,Mn_LPrime_Fit(:),'-x','MarkerSize', 14, 'LineWidth',2)
hold on
semilogy(momentOrder,Mn_LPrime_Predict(:),'-+','MarkerSize', 14, 'LineWidth',2)
semilogy(momentOrder,iterM_L(end,:), '--o','MarkerSize', 15, 'LineWidth',2)
semilogy(momentOrder,iterM_L(floor(ntPrime/2),:), '--o','MarkerSize', 10, 'LineWidth',2)
xlabel('$k$','Interpreter','Latex')
ylabel('$\hat{M}_{k,\mathrm{L}}^{\prime}$','Interpreter','Latex')
legend('$\hat{M}_{k,\mathrm{L}}^{\prime}$ (short-screw)', '$\hat{M}_{k,\mathrm{L}}^{\prime}$ (long-screw)','$\hat{M}_{k,\mathrm{L}}^{\prime}(t^{\prime}=\mathit{MRT}_1)$','$\hat{M}_{k,\mathrm{L}}^{\prime}(t^{\prime}=\mathit{MRT}_2)$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes16,'FontSize', 16)