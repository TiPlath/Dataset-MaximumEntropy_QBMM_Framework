%% Function using DQMOM for twin-screw wet granulation dataset analysis of different distributions
function [steadyStateTime,realSimulationTime,d32,Omega_t,Omega_V] = MRTd32Match(ExperimentalDistribution,a_t,betaPrime_L,MRT)
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define important variables
% The total volume [m^3] of a typical sample from twin-screw wet granulation
V_P = 2.05e-6;
% Total number of particles (zeroth moment)
N_P = 8.78888e+07;
% time step
dtPrime=0.01;
% Number of dirac-delta distributed classes (weights and nodes, 1,...,7)
N_delta=2;
%% read custom Distribution
% call function to read distribution (v_L CDF)
CDFv_L = ReadData('Data/InitialLactoseMCCPVP_CDF-v_L.csv',1e6,0,0.995);
CDFv_L_Fit = ReadData("Data/" + ExperimentalDistribution + ".csv",1e6,0,1);
% Normalize and scale by V_P to represent volume of a typical wet
% granulation sample (assume volume should stay constant)
CDFv_L(:,2) = CDFv_L(:,2)/max(CDFv_L(:,2));
CDFv_L(:,2) = CDFv_L(:,2)*V_P;
CDFv_L_Fit(:,2) = CDFv_L_Fit(:,2)/max(CDFv_L_Fit(:,2));
CDFv_L_Fit(:,2) = CDFv_L_Fit(:,2)*V_P;

%% Conversions into different representations of the distribution function
% convert CDFv_L to PDF n_V for the DQMOM
% Convert CDF v_L to volume-based VDF (v_V)
kappa = pi/6;
CDFv_V(:,2) = CDFv_L(:,2);
CDFv_V(:,1) = kappa .* CDFv_L(:,1).^3;
CDFv_V_Fit(:,2) = CDFv_L_Fit(:,2);
CDFv_V_Fit(:,1) = kappa .* CDFv_L_Fit(:,1).^3;
% Convert CDF v_V to PDF v_V
v_V = convertCDFtoPDF(CDFv_V);
v_V_Fit = convertCDFtoPDF(CDFv_V_Fit);
% Convert volume-based VDF (v_V) to volume-based NDF (n_V)
n_V(:,1) = v_V(:,1);
n_V(:,2) = v_V(:,2)./v_V(:,1);
n_V_Fit(:,1) = v_V_Fit(:,1);
n_V_Fit(:,2) = v_V_Fit(:,2)./v_V_Fit(:,1);

% Conversion to get n_L
% convert PDF n_V to CDF n_V
CDFn_V = convertPDFtoCDF(n_V);
CDFn_V_Fit = convertPDFtoCDF(n_V_Fit);
% convert CDFn_V to CDFn_L
CDFn_L = CDFn_V;
CDFn_L(:,1) = (CDFn_V(:,1)/kappa).^(1/3);
CDFn_L_Fit = CDFn_V_Fit;
CDFn_L_Fit(:,1) = (CDFn_V_Fit(:,1)/kappa).^(1/3);
% convert CDFn_L to PDF n_L
n_L = convertCDFtoPDF(CDFn_L);
n_L_Fit = convertCDFtoPDF(CDFn_L_Fit);

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
% Aggregation rate
% a_t = 5e4;
Omega_a = a_t*Omega_V;
aPrime(1:N_delta,1) = 1;
% Time scale, where aggregation is assumed to be dominant [s]
Omega_t = Omega_P^(-1)*Omega_a^(-1);
% Growth rate length-based
gPrime_L(1:N_delta,1) = 0.00152;
% Breakage rate length-based
betaPrime_L(1:N_delta,1) = betaPrime_L;
% Number of fragments per broken particle length-based
N_fL = 8;
% symmetric fragentation for a length-based PSD (here: n_L)
b_alpha_L = @(L_alpha,k,N_fL) N_fL^((3-k)/3)*L_alpha.^k;

%% determine maximum simulation time and initialise parameters
% max time is similar to 10x MRT
tmaxPrime=10*MRT*Omega_t^(-1);
% Number of time steps
ntPrime=round(tmaxPrime/dtPrime)+1;
% Maximum number of moments
mMax = 2*N_delta;
% Initialize different iteration vectors (useful for debugging and plotting)
% A matrix storing the volume-based momenta of every time step
iterM_L=zeros(ntPrime,2*N_delta);
% store weights and nodes of every iteration
iterL = zeros(ntPrime,N_delta);
iterw_L = zeros(ntPrime,N_delta);

%% Define scaled variables (pre-processing)
% define dimensionless particle Volume length [m] and distribution function n_L [p/m]
% define dimensionless n_L
n_LPrime(:,1) = n_L(:,1)./(Omega_V^(1/3));
n_LPrime(:,2) = n_L(:,2).*Omega_V^(1/3).*Omega_P^(-1);
n_LPrime_Fit(:,1) = n_L_Fit(:,1)./(Omega_V^(1/3));
n_LPrime_Fit(:,2) = n_L_Fit(:,2).*Omega_V^(1/3).*Omega_P^(-1);
% dimensionless time
tPrime=(0:ntPrime-1);
%% evolve abscissae and xweights in time according to DQMOM
% Compute length based moments of our distribution (n_L)
Mn_LPrime = ComputeMoments(n_LPrime(:,1), n_LPrime(:,2), mMax);
Mn_LPrime_Fit = ComputeMoments(n_LPrime_Fit(:,1), n_LPrime_Fit(:,2), mMax);
% we compute weights an nodes of n_L using the Wheeler algorithm
[Ln_LPrime,wn_LPrime] = Wheeler(Mn_LPrime(1:2*N_delta),N_delta);
% Store inital moments, nodes and weights of the Volume based distribution for
% further processing (n_L)
iterM_L(1,:) = getMomenta(Ln_LPrime,wn_LPrime);
% Store inital nodes and weights in iteration vectors
iterL(1,1:N_delta) = Ln_LPrime;
iterw_L(1,1:N_delta) = wn_LPrime;
% Time loop, forward Euler
for i=2:ntPrime
    % Compute one time step of DQMOM (with n_L)
    [Ln_LPrime,wn_LPrime] = DirectQuadratureMethodOfMomentsNonDimensionalizedLengthBased(iterL(i-1,1:N_delta)',iterw_L(i-1,1:N_delta)',dtPrime,gPrime_L,aPrime,betaPrime_L,b_alpha_L,N_fL);
    % compute momenta from weights and nodes (with n_L)
    iterM_L(i,:) = getMomenta(Ln_LPrime,wn_LPrime);

    % Assign weights and nodes to iteration vectors
    iterL(i,1:N_delta) = Ln_LPrime;
    iterw_L(i,1:N_delta) = wn_LPrime;
end

%% Postprocessing

% Find the time at which the function reaches a steady state
tolerance = 1e-4;  % Adjust this tolerance as needed
for i = 2:length(tPrime)
    % Calculate the rate of change
    rateOfChange = abs(iterM_L(i,1) - iterM_L(i - 1,1)) / (tPrime(i) - tPrime(i - 1));
    
    if rateOfChange < tolerance
        steadyStateTime = tPrime(i)*dtPrime;
        break;
    end
end
% precision in percentage
precision = 0.1;
realSimulationTime = steadyStateTime*Omega_t;
if steadyStateTime > 0
    fprintf('Steady state is reached at t'' = %.2f; this is equal to t = %.2f\n', steadyStateTime, realSimulationTime);
    if abs(MRT - realSimulationTime)/MRT * 100 < precision
        % do nothing
    elseif realSimulationTime < MRT
        fprintf(2, 'Decrease aggregation rate to get a higher steady state time (precision %.1f)\n', precision)
    elseif realSimulationTime > MRT
        fprintf(2, 'Increase aggregation rate to get a lower steady state time (precision %.1f)\n', precision)
    end
else
    fprintf('Steady state not reached within the specified tolerance.\n');
end

%% plot d32 fit
stepsize = 10;
d32_N13 = (Mn_LPrime_Fit(4))/(Mn_LPrime_Fit(3)) * Omega_V^(1/3);
d32 = zeros(1,ntPrime);
d32 = (iterM_L(:,4))./(iterM_L(:,3)) * Omega_V^(1/3);

d32Error = (d32(end)-d32_N13)/d32_N13 * 100;
if abs(d32Error) < precision
    fprintf("Sauter mean diameter is fit with an error of: %.2f%%\n\n", d32Error)
elseif d32Error > 0
    fprintf("Sauter mean diameter is fit with an error of: %.2f%%\n", d32Error)
    fprintf(2, 'Increase breakage to get a better fit (precision %.1f)\n\n', precision)
elseif d32Error < 0
    fprintf("Sauter mean diameter is fit with an error of: %.2f%%\n", d32Error)
    fprintf(2, 'Decrease breakage to get a better fit (precision %.1f)\n\n', precision)
end
d32 = d32(end);
end