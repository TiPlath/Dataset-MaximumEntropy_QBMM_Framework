%% Non-dimensional DQMOM solution for special cases (section 4.4 to 4.6)
%
% We compute the non-dimensional DQMOM solution for pure aggregation, pure
% breakage and pure growth using constant kernels. Analytical solutions are
% compared to the DQMOM solution for all available moments.
% We use a mixture of lactose, polyvinylpyrrolidon and microcrystalline cellulose
% from an open-source twin-screw wet granulation dataset (Plath et al.
% 2021)
% The user can choose which special case to run by command line input.

clear all
close all
clc
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% Ask user which rate kernel to compute
computeAggregation = 0;
computeBreakage = 0;
computeGrowth = 0;
reply = input("which rate kernel do you want to compute? \n1 = Aggregation (default)\n2 = Breakage\n3 = Growth\n\n","s");
switch reply
    case "1"
        computeAggregation = 1;
    case "2"
        computeBreakage = 1;
    case "3"
        computeGrowth = 1;  
    otherwise
        computeAggregation = 1;
end
%% Typical Kernel values for non-dimensionalization
% In this script only singlue rate kernels are active, since the
% non-dimensionalisation is based on each individual rate mechanism.
if computeGrowth
    % Growth rate [m^3/s]
    g_0 = 1e-15;
    % max time
    tmaxPrime = 14;
end
if computeAggregation
    % Aggregation rate [1/s]
    a_0 = 1e-8;
    % max time
    tmaxPrime = 200;
end
if computeBreakage
    % Breakage rate[1/s]
    beta_0 = 2e-3;
    % max time
    tmaxPrime = 4;
end
%% define important variables
% The total volume [m^3] of a typical sample from twin-screw wet granulation
V_P = 2.05e-6;
% Total number of particles (zeroth moment)
N_P = 8.78888e+07;
% number of times to plot the MER
ntimes = 4;
% resolution of the reconstructed PSD by MER
resolution = 300;
% time step
dtPrime=0.1;
% symmetric fragmentation for a volume-based PSD (here: n_V)
b_alpha = @(V_alpha,k,N_f) N_f^(1-k)*V_alpha.^k;
% Number of fragments per broken particle
N_f = 2;
% Number of time steps
ntPrime=round(tmaxPrime/dtPrime)+1;
% Number of dirac-delta distributed classes (weights and nodes, 1,...,7)
N_delta=3;
% Maximum number of moments
mMax = 2*N_delta;
% Initialize different iteration vectors (useful for debugging and plotting)
% A matrix storing the volume-based momenta of every time step
iterM_V=zeros(ntPrime,2*N_delta);
%Initialize k-vector to store the realized moments of each MER
iterk = zeros(ntPrime,1);
% store weights and nodes of every iteration
iterV = zeros(ntPrime,N_delta);
iterw_V = zeros(ntPrime,N_delta);
%% read custom Distribution
% call function to read distribution (v_L CDF)
CDFv_L = ReadData('Data/InitialLactoseMCCPVP_CDF-v_L.csv',1e6,0,0.995);
% Normalize and scale by V_P to represent volume of a typical wet
% granulation sample
CDFv_L(:,2) = CDFv_L(:,2)/max(CDFv_L(:,2));
CDFv_L(:,2) = CDFv_L(:,2)*V_P;
% Initialize x-vector to store domain data of every iteration
iterx = zeros(ntPrime,resolution);
% store PSD data of every iteration
iterPSD = zeros(ntPrime,resolution);
%% Density distribution conversions
% Convert CDF v_L to PDF v_L for the MER
v_L = convertCDFtoPDF(CDFv_L);
% Convert CDF v_L to volume-based CDF (v_V)
k_v = pi/6;
CDFv_V(:,2) = CDFv_L(:,2);
CDFv_V(:,1) = k_v .* CDFv_L(:,1).^3;
% Convert CDF v_V to PDF v_V
v_V = convertCDFtoPDF(CDFv_V);
% Convert volume-based VDF (v_V) to volume-based NDF (n_V)
n_V(:,1) = v_V(:,1);
n_V(:,2) = v_V(:,2)./v_V(:,1);

%% Define Scale variables for non-dimensionalisation at t=0 (pre-processing)
% Compute total number of particles to define particle unit scale P
M_0 =  sum(n_V(1:end-1,2).*diff(n_V(:,1)));
% Compute total volume to define length scale L
M_1V = sum(n_V(1:end-1,1).*n_V(1:end-1,2).*diff(n_V(:,1)));
% Particle scale [p]
Omega_P = M_0;
M0Prime = 1;
% Volume scale [m^3]
Omega_V = (M_1V * Omega_P^(-1));
M1Prime = 1;
% Time scale depending on rate mechanism [s]
if computeGrowth ~= 0
    Omega_t = g_0^(-1) * Omega_V;
    gPrime(1:N_delta,1) = 1;
    aPrime(1:N_delta,1) = 0;
    betaPrime(1:N_delta,1) = 0;
end
if computeAggregation ~= 0
    Omega_t = a_0^(-1) * Omega_P^(-1);
    gPrime(1:N_delta,1) = 0;
    aPrime(1:N_delta,1) = 1;
    betaPrime(1:N_delta,1) = 0;
end
if computeBreakage ~= 0
    Omega_t = beta_0^(-1);
    gPrime(1:N_delta,1) = 0;
    aPrime(1:N_delta,1) = 0;
    betaPrime(1:N_delta,1) = 1;
end
if (computeBreakage(1) == 0 && computeAggregation(1) == 0 && computeGrowth(1) == 0)
    Omega_t = 1;
    gPrime(1:N_delta,1) = 0;
    aPrime(1:N_delta,1) = 0;
    betaPrime(1:N_delta,1) = 0;
end

%% Define scaled variables (pre-processing)
% define dimensionless particle Volume V [m^3] and distribution function n_V [p/m^3]
% to ensure we get dimensionless moments M_V
n_VPrime(:,1) = n_V(:,1)./(Omega_V);
n_VPrime(:,2) = n_V(:,2).*Omega_V.*Omega_P^(-1);
% define dimensionless particle Length L [m] and distribution function v_L [p m^3/m]
% to ensure we get dimensionless moments Mv_L
v_LPrime(:,1) = v_L(:,1)./(Omega_V)^(1/3);
v_LPrime(:,2) = v_L(:,2).*Omega_V^(-2/3).*Omega_P^(-1);
% dimensionless time
tPrime=(0:ntPrime-1);
%% get initial weights and nodes for a volume-based distribution (v_L)
n_LPrime(:,1) = v_LPrime(:,1);
n_LPrime(:,2) = v_LPrime(:,2)./(k_v.*v_LPrime(:,1).^3);
% Compute length-based moments of distribution for the initial reconstruction
% v_L
Mv_LPrime = ComputeMoments(v_LPrime(:,1), v_LPrime(:,2), mMax);
% we compute weights and nodes using the Wheeler algorithm
% v_L
[Lv_LPrime,wv_LPrime] = Wheeler(Mv_LPrime(1:2*N_delta),N_delta);
% reconstruct the initial density distribution by maximum entropy approach
[d,E_V,M_,initx,initPSD,c,lambda] = plotInitLeastErrorPSD(v_LPrime(:,1),v_LPrime(:,2),Lv_LPrime,Mv_LPrime,mMax);
% fix xMin to the minimum particle size of the initial reconstruction
xMin = 0;
% store initial densities, reconstruction domain to the iteration vectors
iterx(1,1:end) = initx;
iterPSD(1,1:end) = initPSD;
% Legend for plotting
str = {'$v_{\mathrm{L}}^{\prime}(t^{\prime} = 0)$', '$\bar{v}_{\mathrm{L}}^{\prime}(t^{\prime} = 0)$'};
%% evolve nodes and weights in time according to DQMOM
% Compute volume-based moments of distribution (n_V)
M_VPrime = ComputeMoments(n_VPrime(:,1), n_VPrime(:,2), mMax);
% we compute weights and nodes of n_V using the Wheeler algorithm 
[VPrime,w_VPrime] = Wheeler(M_VPrime(1:2*N_delta),N_delta);

% Store inital moments, nodes and weights of the Volume based distribution for
% further processing (n_V)
iterM_V(1,:) = getMomenta(VPrime,w_VPrime);
% Store inital nodes and weights in iteration vectors
iterV(1,1:N_delta) = VPrime;
iterw_V(1,1:N_delta) = w_VPrime;
% Time loop, forward Euler
for i=2:ntPrime
    % Compute one time step of DQMOM (with n_V)
    [VPrime,w_VPrime] = DirectQuadratureMethodOfMomentsNonDimensionalized(iterV(i-1,1:N_delta)',iterw_V(i-1,1:N_delta)',dtPrime,gPrime,aPrime,betaPrime,b_alpha,N_f);
    % compute momenta from weights and nodes (with n_V)
    iterM_V(i,:) = getMomenta(VPrime,w_VPrime);
    % Assign weights and nodes to iteration vectors
    iterV(i,1:N_delta) = VPrime;
    iterw_V(i,1:N_delta) = w_VPrime;
    if mod(i-1,round(ntPrime/ntimes))==0
        % convert volume-based number nodes and weights (n_V) to length-based
        % volume nodes and weights (v_L)
        Lv_LPrime = (iterV(i,1:N_delta)./k_v).^(1/3);
        wv_LPrime = iterw_V(i,1:N_delta).*iterV(i,1:N_delta);
        % set fixed distance to extend over the actual function
        % (c*L_{N_{\delta}})
        dxMax = abs(max(Lv_LPrime)*(c-1));
        % Reconstruct and plot the PSD from the distributions moments
        [PSD,x,k,lambda] = plotMaximumEntropyReconstruction(Lv_LPrime,wv_LPrime,mMax,xMin,dxMax,v_LPrime(:,2));
        xlabel('$L^{\prime}$','Interpreter','Latex')
        ylabel('$v_{\mathrm{L}}^{\prime}$','Interpreter','Latex')        
        str = [str , sprintf('$\\bar{v}_{\\mathrm{L}}^{\\prime} (t^{\\prime} = %.2f)$' , tPrime(i)*dtPrime)];
        legend(str, 'location', 'best','Interpreter','Latex')
        % Store densities, domain and number of moments used in the
        % iteration vectors
        iterPSD(i,1:length(PSD)) = PSD;
        iterx(i,1:length(x)) = x;
        iterk(i) = k;
    end
end
%% Post-Processing
%% Define analytical moment solutions for constant growth, aggregation and Breakage
% Growth
if computeGrowth ~= 0
    syms M_0Growth(t_)
    M_0Growth(t_) = iterM_V(1,1);
    M_1Growth(t_) = iterM_V(1,2) + gPrime(1,1)*iterM_V(1,1)*t_;
    M_2Growth(t_) = iterM_V(1,3) + 2*gPrime(1,1)*iterM_V(1,2)*t_ + gPrime(1,1)^2*iterM_V(1,1)*t_.^2;
    M_3Growth(t_) = iterM_V(1,4) + 3*gPrime(1,1)*iterM_V(1,3)*t_ + 3*gPrime(1,1)^2*iterM_V(1,2)*t_.^2 + gPrime(1,1)^3*iterM_V(1,1)*t_.^3;
    M_4Growth(t_) = iterM_V(1,5) + 4*gPrime(1,1)*iterM_V(1,4)*t_ + 6*gPrime(1,1)^2*iterM_V(1,3)*t_.^2 + gPrime(1,1)^3*4*iterM_V(1,2)*t_.^3 + gPrime(1,1)^4*iterM_V(1,1)*t_.^4;
    M_5Growth(t_) = iterM_V(1,6) + 5*gPrime(1,1)*iterM_V(1,5)*t_ + 10*gPrime(1,1)^2*iterM_V(1,4)*t_.^2 + 2*gPrime(1,1)^3*5*iterM_V(1,3)*t_.^3 + gPrime(1,1)^4*5*iterM_V(1,2)*t_.^4 + gPrime(1,1)^5*iterM_V(1,1)*t_.^5;
end

% Aggregation
if computeAggregation ~=0
    syms M_0Agg(t_) M_1Agg(t_)
    M_0Agg(t_) = 2*iterM_V(1,1)/(2+aPrime(1,1)*iterM_V(1,1)*t_);
    M_1Agg(t_) = iterM_V(1,2);
    M_2Agg(t_) = iterM_V(1,3)+ aPrime(1,1)*iterM_V(1,2)^2*t_;
    M_3Agg(t_) = iterM_V(1,4)+ 3*aPrime(1,1)*iterM_V(1,2)*iterM_V(1,3)*t_ + 3/2*aPrime(1,1)^2*iterM_V(1,2)^3*t_.^2;
    M_4Agg(t_) = iterM_V(1,5)+ 3*aPrime(1,1)*iterM_V(1,3)^2*t_ + 4*aPrime(1,1)*iterM_V(1,2)*iterM_V(1,4)*t_ + 9*aPrime(1,1)^2*iterM_V(1,2)^2*iterM_V(1,3)*t_.^2 + 3*aPrime(1,1)^3*iterM_V(1,2)^4*t_.^3;
    M_5Agg(t_) = iterM_V(1,6)+ 10*aPrime(1,1)*iterM_V(1,3)*iterM_V(1,4)*t_ + 5*aPrime(1,1)*iterM_V(1,2)*iterM_V(1,5)*t_ + 45/2*aPrime(1,1)^2*iterM_V(1,3)^2*iterM_V(1,2)*t_.^2 + 15*aPrime(1,1)^2*iterM_V(1,2)^2*iterM_V(1,4)*t_.^2 + 30*aPrime(1,1)^3*iterM_V(1,2)^3*iterM_V(1,3)*t_.^3 + 15/2*aPrime(1,1)^4*iterM_V(1,2)^5*t_.^4;
end 

% Breakage
if computeBreakage ~=0
    M_kBreakage = @(k,t) exp((2^(1-k)-1)*betaPrime(1,1).*t) .* iterM_V(1,k+1);
end

if (computeBreakage(1) == 0 && computeAggregation(1) == 0 && computeGrowth(1) == 0)
    for i = 1:mMax
        M(i) = iterM_V(1,i);
    end
end
%% plot solutions
stepsize = round(ntPrime*dtPrime);
% compute error E_V
% aggregation
if computeAggregation~=0
   for i = 1:length(iterM_V)
        E_V(i,1) = abs(iterM_V(i,1) - M_0Agg((i-1)*dtPrime))/M_0Agg((i-1)*dtPrime);
        E_V(i,2) = abs(iterM_V(i,2) - M_1Agg((i-1)*dtPrime))/M_1Agg((i-1)*dtPrime);
        E_V(i,3) = abs(iterM_V(i,3) - M_2Agg((i-1)*dtPrime))/M_2Agg((i-1)*dtPrime);
        E_V(i,4) = abs(iterM_V(i,4) - M_3Agg((i-1)*dtPrime))/M_3Agg((i-1)*dtPrime);
        E_V(i,5) = abs(iterM_V(i,5) - M_4Agg((i-1)*dtPrime))/M_4Agg((i-1)*dtPrime);
        E_V(i,6) = abs(iterM_V(i,6) - M_5Agg((i-1)*dtPrime))/M_5Agg((i-1)*dtPrime);
    end
end
% growth
if computeGrowth~=0
    for i = 1:length(iterM_V)
        E_V(i,1) = abs(iterM_V(i,1) - M_0Growth((i-1)*dtPrime))/M_0Growth((i-1)*dtPrime);
        E_V(i,2) = abs(iterM_V(i,2) - M_1Growth((i-1)*dtPrime))/M_1Growth((i-1)*dtPrime);
        E_V(i,3) = abs(iterM_V(i,3) - M_2Growth((i-1)*dtPrime))/M_2Growth((i-1)*dtPrime);
        E_V(i,4) = abs(iterM_V(i,4) - M_3Growth((i-1)*dtPrime))/M_3Growth((i-1)*dtPrime);
        E_V(i,5) = abs(iterM_V(i,5) - M_4Growth((i-1)*dtPrime))/M_4Growth((i-1)*dtPrime);
        E_V(i,6) = abs(iterM_V(i,6) - M_5Growth((i-1)*dtPrime))/M_5Growth((i-1)*dtPrime);
    end
end
% breakage
if computeBreakage~=0
    for j = 0:mMax-1
        for i = 1:length(iterM_V)
            E_V(i,j+1) = abs(iterM_V(i,j+1) - M_kBreakage(j,(i-1)*dtPrime))/M_kBreakage(j,(i-1)*dtPrime);
        end
    end
end

% none
if (computeBreakage == 0 && computeAggregation == 0 && computeGrowth == 0)
    for j = 0:mMax-1
        for i = 1:length(iterM_V)
            E_V(i,j+1) = abs(iterM_V(i,j+1) - M(j+1))/M(j+1);
        end
    end
end

% plot error between DQMOM and analytical solution
figure(5)
axes5 = axes('Parent', figure(5));
for i = 1:size(E_V,2)
    plot(tPrime(1:stepsize:length(E_V))*dtPrime,E_V(1:stepsize:end,i), "-",'MarkerSize', 14, 'LineWidth',2)
    hold on
end
xlabel('$t^{\prime}$','Interpreter','Latex')
ylabel('$E$ [\%]','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes5,'FontSize', 16)
%% plot DQMOM moments against analytical moments
figure(7)
axes7 = axes('Parent', figure(7));
for i = 1:size(iterM_V,2)
    plot(tPrime(1:stepsize:ntPrime+1)*dtPrime,iterM_V(1:stepsize:end,i)/max(iterM_V(:,i)),'x', 'MarkerSize', 14,  'LineWidth', 2)
    hold on
end
% Reset the color order to repeat the defined colors
set(axes7, 'ColorOrderIndex', 1);
% growth 
if computeGrowth ~= 0
    fplot(subs(M_0Growth)/max(subs(M_0Growth([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2);
    hold on
    fplot(subs(M_1Growth)/max(subs(M_1Growth([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_2Growth)/max(subs(M_2Growth([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_3Growth)/max(subs(M_3Growth([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_4Growth)/max(subs(M_4Growth([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_5Growth)/max(subs(M_5Growth([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    legend("$\widetilde{M}_{0}^{g\prime}$","$\widetilde{M}_{1,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{2,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{3,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{4,\mathrm{V}}^{g\prime}$","$\widetilde{M}_{5,\mathrm{V}}^{g\prime}$",'location','best', 'FontSize', 14,'Interpreter','Latex')
end
% aggregation
if computeAggregation ~= 0
    fplot(subs(M_0Agg)/max(subs(M_0Agg([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2);
    hold on
    fplot(subs(M_1Agg)/max(subs(M_1Agg([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_2Agg)/max(subs(M_2Agg([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_3Agg)/max(subs(M_3Agg([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_4Agg)/max(subs(M_4Agg([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    fplot(subs(M_5Agg)/max(subs(M_5Agg([0,tPrime(1:ntPrime)*dtPrime]))),[0 tmaxPrime], 'LineWidth', 2)
    legend("$\widetilde{M}_{0}^{a\prime}$","$\widetilde{M}_{1,\mathrm{V}}^{a\prime}$","$\widetilde{M}_{2,\mathrm{V}}^{a\prime}$","$\widetilde{M}_{3,\mathrm{V}}^{a\prime}$","$\widetilde{M}_{4,\mathrm{V}}^{a\prime}$","$\widetilde{M}_{5,\mathrm{V}}^{a\prime}$","$M_{k,\mathrm{V}}^{a\prime}$",'location','best', 'FontSize', 14,'Interpreter','Latex')
end
% breakage
if computeBreakage ~= 0
    for i = 0:mMax-1
        plot(tPrime(1:ntPrime)*dtPrime,M_kBreakage(i,tPrime(1:ntPrime)*dtPrime)/max(M_kBreakage(i,tPrime(1:ntPrime)*dtPrime)), 'LineWidth', 2)
        hold on
    end
    legend("$\widetilde{M}_{0}^{b \prime}$","$\widetilde{M}_{1,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{2,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{3,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{4,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{5,\mathrm{V}}^{b\prime}$","$M_{k,\mathrm{V}}^{b\prime}$",'location','best', 'FontSize', 14,'Interpreter','Latex')
end
% none
if (computeBreakage == 0 && computeAggregation == 0 && computeGrowth == 0)
   for i = 1:mMax
        plot(tPrime(1:ntPrime)*dtPrime,(M(i)*tPrime)./(M(i)*tPrime), 'LineWidth', 2)
        hold on
    end
    legend("$\widetilde{M}_{0}^{b \prime}$","$\widetilde{M}_{1,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{2,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{3,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{4,\mathrm{V}}^{b\prime}$","$\widetilde{M}_{5,\mathrm{V}}^{b\prime}$","$M_{k,\mathrm{V}}^{b\prime}$",'location','best', 'FontSize', 14,'Interpreter','Latex')
end

xlabel('$t^\prime$','Interpreter','Latex')
ylabel('$\hat{M}_{k,\mathrm{V}}^\prime$','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes7,'FontSize', 16)