%% Application to twin-screw wet granulation data: Data extrapolation (section 5.3)

% We extrapolate our created model to datapoints outside the range of the 
% twin-screw wet granulation dataset (Plath et al. 2021) with the lactose,
% polyvinylpyrrolidon and microcrystalline cellulose mixture.
% 1. For a certain SFL and L/S parameter pair, we predict the mean residence 
%    time MRT as well as the functions found in the dataset analysis (1) to
%    predict the Sauter mean diameter d_32, minimum aggregation rate a_s 
%    and non-dimensional breakage rate beta_t'.
% 2. We use the same kernel values to predict the long-screw process by 
%    simulating for another mean residence time.

clc
close all
clear variables
% disable legend warnings
warningState = warning('off', 'MATLAB:legend:IgnoringExtraEntries');
warning('off', 'MATLAB:nearlySingularMatrix');
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define important variables (pre-processing)
% The total volume [m^3] of a typical sample from twin-screw wet granulation
V_P = 2.05e-6;
% Total number of particles (zeroth moment)
N_P = 8.78888e+07;
% Density
rho_P = 1.540;
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
% LOW PARAMETERS
SFL = 0.02;
LS = 0.05;
% HIGH PARAMETERS
% SFL = 0.2;
% LS = 0.4;
% Number of fragments per broken particle length-based
N_fL =8;
%% read custom Distribution
% call function to read distribution (v_L CDF)
CDFv_L = ReadData('Data/InitialLactoseMCCPVP_CDF-v_L.csv',1e6,0,0.995);
% Normalize and scale by V_P to represent volume of a typical wet
% granulation sample (assume volume should stay constant)
CDFv_L(:,2) = CDFv_L(:,2)/max(CDFv_L(:,2));
CDFv_L(:,2) = CDFv_L(:,2)*V_P;
%% Conversions into different representations of the distribution function
% Convert CDF v_L to PDF v_L for the MER
v_L = convertCDFtoPDF(CDFv_L);
% convert CDFv_L to PDF n_V for the DQMOM
% Convert CDF v_L to volume-based VDF (v_V)
kappa = pi/6;
CDFv_V(:,2) = CDFv_L(:,2);
CDFv_V(:,1) = kappa .* CDFv_L(:,1).^3;
% Convert CDF v_V to PDF v_V
v_V = convertCDFtoPDF(CDFv_V);
% Convert volume-based VDF (v_V) to volume-based NDF (n_V)
n_V(:,1) = v_V(:,1);
n_V(:,2) = v_V(:,2)./v_V(:,1);
% Conversion to get n_L
% convert PDF n_V to CDF n_V
CDFn_V = convertPDFtoCDF(n_V);
% convert CDFn_V to CDFn_L
CDFn_L = CDFn_V;
CDFn_L(:,1) = (CDFn_V(:,1)/kappa).^(1/3);
% convert CDFn_L to PDF n_L
n_L = convertCDFtoPDF(CDFn_L);
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

%% predictions for MRT, d_32, a_s and beta_tPrime from dataset analysis
% mean residence time prediction from SFL and LS short-screw
MRTFunc = @(params, xy) params(1) * xy(:,1)+ params(2) * xy(:,2) + params(3);
abc1 = [123.46,45.16,3.89];
parameterPair = [SFL, LS];
MRT = MRTFunc(abc1,parameterPair);

% Sauter mean diameter prediction from SFL and LS
d32Func = @(params,xy) (params(1).*xy(:,2).^2) .* (xy(:,1)./rho_P).^(1/3);
abc2 = 0.0605;
d32 = d32Func(abc2,parameterPair);

% minimum aggregation rate prediction from MRT, d32
t_opt = 7.5;
a_sFunc = @(params, xy) params(1) * ((xy(:,1)-t_opt)/t_opt).^(params(2)) .* (xy(:,2)*Omega_V^(-1/3)).^(params(3));
% a_c, p_t, p_m
abc3 = [0.3184, -0.55, -0.181];
a_s = a_sFunc(abc3, [MRT,d32]);  
a_t = a_s*(Omega_P*Omega_V)^(-1);

% non-dimensional breakage rate from d32#
beta_tPrimeFunc = @(params, x) params(1).*x.^(-params(2));
abc4 = [2.625, 0.914];
d32Prime = d32*Omega_V^(-1/3);
beta_tPrime = beta_tPrimeFunc(abc4,d32Prime);
%% set rates for the kernels using the predictions and compute timescale
% Aggregation rate
% a_t = 5e4;
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
% max time is similar to 10x MRT
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

%% Define scaled variables (pre-processing)
% define dimensionless particle Volume V [m^3] and distribution function n_V [p/m^3]
% to ensure we get dimensionless moments M_V
n_VPrime(:,1) = n_V(:,1)./(Omega_V);
n_VPrime(:,2) = n_V(:,2).*Omega_V.*Omega_P^(-1);
% define dimensionless particle Length L [m] and distribution function v_L [p m^3/m]
% to ensure we get dimensionless moments Mv_L
v_LPrime(:,1) = v_L(:,1)./(Omega_V)^(1/3);
v_LPrime(:,2) = v_L(:,2).*Omega_V^(-2/3).*Omega_P^(-1);
% define dimensionless n_L
n_LPrime(:,1) = n_L(:,1)./(Omega_V^(1/3));
n_LPrime(:,2) = n_L(:,2).*Omega_V^(1/3).*Omega_P^(-1);
% dimensionless time
tPrime=(0:ntPrime-1)*dtPrime;

%% get initial weights and nodes for a volume-based distribution (v_L)
% Compute length-based moments of distribution for the initial reconstruction
Mv_LPrime = ComputeMoments(v_LPrime(:,1), v_LPrime(:,2), mMax);
% we compute weights and nodes using the Wheeler algorithm
[Lv_LPrime,wv_LPrime] = Wheeler(Mv_LPrime(1:2*N_delta),N_delta);
% plot the initial delta-pdf by ME approach using moments of v_L
[d,initE,M_,initx,initPSD,c,lambda] = plotInitLeastErrorPSD(v_LPrime(:,1),v_LPrime(:,2),Lv_LPrime,Mv_LPrime,mMax);
% plot final distributions
figure(1)
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlim([0 1.01e3])
% fix xMin to the minimum particle size of the initial reconstruction
xMin = 0;
% store initial densities, reconstruction domain to the iteration vectors
iterx(1,1:end) = initx;
iterPSD(1,1:end) = initPSD;
iterLambda(1,1:length(lambda)) = lambda;
%% evolve abscissae and weights in time according to DQMOM
% n_L
Mn_LPrime = ComputeMoments(n_LPrime(:,1), n_LPrime(:,2), mMax);
% We compute weights and nodes from n_L using the Wheeler algorithm
[Ln_LPrime,wn_LPrime] = Wheeler(Mn_LPrime(1:2*N_delta),N_delta);
% Store inital moments, nodes and weights of the Volume based distribution for
% further processing (n_L)
iterM_L(1,:) = getMomenta(Ln_LPrime,wn_LPrime);
% Store inital nodes and weights in iteration vectors
iterL(1,1:N_delta) = Ln_LPrime;
iterw_L(1,1:N_delta) = wn_LPrime;
iterdxMax(1) = abs(max(Ln_LPrime)*(c-1));
% define counter
o = 0;
str = {'$v_{\mathrm{L}}(t = 0)$', '$\bar{v}_{\mathrm{L}}(t = 0)$'};
% Time loop, forward Euler
for i=2:ntPrime
    % Compute one time step of DQMOM (n_L)
    [Ln_LPrime,wn_LPrime] = DirectQuadratureMethodOfMomentsNonDimensionalizedLengthBased(iterL(i-1,1:N_delta)',iterw_L(i-1,1:N_delta)',dtPrime,gPrime_L,aPrime,betaPrime_L,b_alpha_L,N_fL);
    % compute momenta from weights and nodes (n_L)
    iterM_L(i,:) = getMomenta(Ln_LPrime,wn_LPrime);
    % Assign weights and nodes to iteration vectors
    iterL(i,1:N_delta) = Ln_LPrime;
    iterw_L(i,1:N_delta) = wn_LPrime;
    % compute distance to extend over the actual function
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
            str = [str , sprintf('$\\bar{v}_{\\mathrm{L}}^{\\prime} (t^{\\prime} = \\mathit{MRT})$' , tPrime(i))];
            iterPSD(2,1:length(PSD)) = PSD;
            iterx(2,1:length(x)) = x;
            iterLambda(2,1:length(lambda)) = lambda;
        end
        if o == 2
            str = [str , sprintf('$\\bar{v}_{\\mathrm{L}}^{\\prime} (t^{\\prime} = 2\\mathit{MRT})$' , tPrime(i))];
            iterPSD(3,1:length(PSD)) = PSD;
            iterx(3,1:length(x)) = x;
            iterLambda(3,1:length(lambda)) = lambda;
        end
        legend(str, 'location', 'best','Interpreter','Latex')
        % Store densities, domain and number of moments used in the,
        % iteration vectors
        iterk(i) = k;
    end    
end

%% Postprocessing
%% scaling the moments to the same units
for k = 1:mMax-1
    iterM_L(:,k+1) = iterM_L(:,k+1).^(1/(k));
end
%% plot moments against moment order semilogarithmic L

momentOrder = 0:1:(size(iterM_L,2) - 1);
figure(16)
axes16 = axes('Parent', figure(16));
semilogy(momentOrder,iterM_L(end,:), '--o','MarkerSize', 15, 'LineWidth',2)
hold on
semilogy(momentOrder,iterM_L(floor(ntPrime/2),:), '--o','MarkerSize', 10, 'LineWidth',2)
xlabel('Moment order','Interpreter','Latex')
ylabel('$\hat{M}_{k_{\mathrm{L}}}^{\prime}$','Interpreter','Latex')
legend('$\hat{M}_{k_{\mathrm{L}}}^{\prime}(t^{\prime}=\mathit{MRT})$','$\hat{M}_{k_{\mathrm{L}}}^{\prime}(t^{\prime}=2\mathit{MRT})$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes16,'FontSize', 16)