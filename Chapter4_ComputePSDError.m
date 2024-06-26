%% Compare the moment error of length-based and volume-based DQMOM (section 4.2)
%
% We compute a volume-based (using n_V) and a length-based (using n_L)
% DQMOM and compare the moment errrors for different timesteps.
% We use a mixture of lactose, polyvinylpyrrolidon and microcrystalline 
% cellulose from an open-source twin-screw wet granulation dataset (Plath 
% et al. 2021)

clear all
close all
clc
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% read custom Distribution
% call function to read distribution (v_L CDF)
CDFv_L = ReadData('Data/InitialLactoseMCCPVP_CDF-v_L.csv',1e6,0,0.995);
% Normalize
CDFv_L(:,2) = CDFv_L(:,2)/100;
%% Define variables (pre-processing)
%%%%%%%%%%%%%%%%%% User defined variables %%%%%%%%%%%%%%%%%%
% resolution of the reconstructed PSD from MER
resolution = 300;
% time step
dt=[1,0.1,0.01];
% max time
tmax=100;
% Number of dirac-delta distributed classes (weights and nodes, 1,...,7)
N_delta=3;
% Growth rate [m^3/s] 
g = ones(N_delta,1);
g(:,1) = 0;
% Aggregation rate [m^3/s]
a = zeros(N_delta,1);
a(:,1) = 1e-15;
% Breakage rate[m^3/s]
beta = zeros(N_delta,1);
beta(:,1) = 0;
% symmetric fragentation for a volume-based PSD (here: n_V)
b_alpha = @(V_alpha,k) 2^(1-k)*V_alpha.^k;
% symmetric fragentation for a length-based PSD (here: n_L)
b_alpha_L = @(L_alpha,k) 2^((3-k)/3)*L_alpha.^k;

figure(5)
axes5 = axes('Parent', figure(5));
hold on
% axis tight
for l = 1:length(dt)
% Number of time steps
nt=round(tmax/dt(l))+1;
% Time
t=(0:nt)*dt(l);
% Maximum number of moments
mMax = 2*N_delta;
%Initialize initial values for iteration vectors (useful for debugging and plotting)
% A matrix storing the volume-based momenta of every time step
iterm_V=zeros(nt,2*N_delta);
% A matrix storing the length-based momenta of every time step
iterm_L=zeros(nt,2*N_delta);
%Initialize k-vector to store the realized moments of each MER
iterk = zeros(nt,1);
% Initialize x-vector to store domain data of every iteration
iterx = zeros(nt,resolution);
% store PSD data of every iteration
iterPSD = zeros(nt,resolution);
% store weights and nodes of every iteration
iterV = zeros(nt,N_delta);
iterw_V = zeros(nt,N_delta);
iterL = zeros(nt,N_delta);
iterw_L = zeros(nt,N_delta);

%% get initial weights and nodes for a volume weighted distribution (v_L)
% Convert CDF v_L to PDF v_L to compute the moments
v_L = convertCDFtoPDF(CDFv_L);
% Compute length-based moments of distribution for the initial reconstruction
Mv_L = ComputeMoments(v_L(:,1), v_L(:,2), mMax);
% we compute weights and nodes using the Wheeler algorithm
[Lv_L,wv_L] = Wheeler(Mv_L(1:2*N_delta),N_delta);
%% Density distribution conversions
% Convert to volume-based VDF (v_V)
k_v = pi/6;
CDFv_V(:,2) = CDFv_L(:,2);
CDFv_V(:,1) = k_v .* CDFv_L(:,1).^3;
% Convert CDF v_V to PDF v_V
v_V = convertCDFtoPDF(CDFv_V);
% Convert volume-based VDF (v_V) to volume-based NDF (n_V)
n_V(:,1) = v_V(:,1);
n_V(:,2) = v_V(:,2)./v_V(:,1);
% Compute volume-based moments of distribution
M_V = ComputeMoments(n_V(:,1), n_V(:,2), mMax);
% we compute weights and nodes using the Wheeler algorithm
[V,w_V] = Wheeler(M_V(1:2*N_delta),N_delta);
% Store inital moments, nodes and weights of the Volume based distribution for
% further processing (n_V)
iterm_V(1,:) = getMomenta(V,w_V);
% Convert n_V to n_L to store length-based moments
% Convert n_V to CDFn_V
CDFn_V = convertPDFtoCDF(n_V);
% convert CDFn_V to CDFn_L
CDFn_L(:,1) = (CDFn_V(:,1)/k_v).^(1/3);
CDFn_L(:,2) = CDFn_V(:,2);
n_L = convertCDFtoPDF(CDFn_L);
% Compute length-based moments of n_L
M_L = ComputeMoments(n_L(:,1), n_L(:,2), mMax);
% Compute weights and nodes of n_L using the Wheeler algorithm
[L,w_L] = Wheeler(M_L(1:2*N_delta),N_delta);
% Store moments of n_L in iteration vector
iterm_L(1,:) = getMomenta(L,w_L);
% Store nodes and weights in iteration vectors
iterV(1,1:N_delta) = V;
iterL(1,1:N_delta) = L;
iterw_V(1,1:N_delta) = w_V;
iterw_L(1,1:N_delta) = w_L;
%% Evolve nodes and weights in time according to DQMOM
% Time loop, forward Euler
for i=2:nt
    % Compute one time step of DQMOM (with n_V)
    [V,w_V] = DirectQuadratureMethodOfMoments(iterV(i-1,1:N_delta)',w_V,dt(l),g,a,beta, b_alpha);
    % Compute one time step of DQMOM (with n_L)
    [L,w_L] = DirectQuadratureMethodOfMomentsLengthBased(iterL(i-1,1:N_delta)',w_L,dt(l),g,a,beta, b_alpha_L);
    % compute momenta from weights and nodes (with n_V)
    iterm_V(i,:) = getMomenta(V,w_V);
    % compute momenta from weights and nodes (with n_L)
    iterm_L(i,:) = getMomenta(L,w_L);
    % Assign weights, nodes, new minimum and maximum particle size to iteration vectors
    iterV(i,1:N_delta) = V;
    iterw_V(i,1:N_delta) = w_V;
    iterL(i,1:N_delta) = L;
    iterw_L(i,1:N_delta) = w_L;
end

%% plot error of Moments
stepsize = round(nt/10);

% compute error
E_L = zeros(nt,1);
% compute E_V only for dt=1
if dt(l)==1
    E_V = zeros(nt,1);

    for i = 1:length(iterm_V)
        E_V(i,1) = abs(iterm_V(i,2) - iterm_V(1,2))/iterm_V(1,2);
    end
    % plot volume-based error
    for i = 1:size(E_V,2)
        plot([0, t(stepsize+1:stepsize:length(E_V))],E_V(1:stepsize:end,1), "x",'MarkerSize', 14, 'LineWidth',2)
    end
end
% compute length-based error for all timesteps
for i = 1:length(iterm_L)
    E_L(i,1) = abs(iterm_L(i,4) - iterm_L(1,4))/iterm_L(1,4);
end
% plot length-based error
for i = 1:size(E_L,2)
    plot([0, t(stepsize+1:stepsize:length(E_L))],E_L(1:stepsize:end,1), "o",'MarkerSize', 14, 'LineWidth',2)
end

end
% options for plotting
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$E_{k,\,\,\,}$ [\%]','Interpreter','Latex')
text(-11.949392265193822,0.000118996305415,-1,char(958),"Interpreter","tex","Rotation",90,"FontSize",13,"FontName","Cambria Math","Color", 0.15*ones(1,3))
legend('$E_{1,\mathrm{V}}, dt=1$','$E_{3,\mathrm{L}}, dt=1$','$E_{3,\mathrm{L}}, dt=0.1$','$E_{3,\mathrm{L}}, dt=0.01$', 'location','best', 'FontSize', 14,'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes5,'FontSize', 16)