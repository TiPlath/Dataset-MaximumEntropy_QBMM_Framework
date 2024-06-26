%% Initial MER reconstruction step (chapter 3)
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% We test the initial maximum entropy reconstruction on a log-normal distribution
% function. A bisection method determines the optimal c-value that minimises 
% the Jenson-Shannon divergence between the reconstruction and the initial 
% distribution of our wet granulation process. The c-value determines the 
% relative distance between the largest node and the maximum particle size 
% of the reconstruction domain.

clear all
close all
clc

%% Read Analytic distribution
% total volume
V_P = 1;
% Log-Normal distribution
f_ = @(x,sigma,mu) (1./((sqrt(2*pi).*sigma.*x))) .* exp((-(log(x)-mu).^2 )./(2*sigma^2));
%% test functions
% Unimodal Log-Normal function
x_f = linspace(1, 150, 150)';
f_uni = V_P * f_(x_f,1/4,3.5);
% Bimodal Log-Normal distribution
f_bi = V_P/2 * f_(x_f,1/3,4) + V_P/2 * f_(x_f,1/7,4.5);

%% Define variables (pre-processing)
% Number of dirac-delta distributed classes (weights and nodes, 1,...,7)
N_uni=3;
N_bi=5;
% Maximum number of moments
mMax_uni = 2*N_uni;
mMax_bi = 2*N_bi;

%% Reconstruct length-based volume distribution (V_L)
% We compute moments of distribution for the initial reconstruction
M_uni = ComputeMoments(x_f, f_uni, mMax_uni);
M_bi = ComputeMoments(x_f, f_bi, mMax_bi);
% we compute weights and nodes using the Wheeler algorithm
[xi_uni,w_uni] = Wheeler(M_uni(1:2*N_uni),N_uni);
[xi_bi,w_bi] = Wheeler(M_bi(1:2*N_bi),N_bi);
% plot the initial delta-pdf by ME approach

% [R_uni,E_uni,M_uni_,initx_uni,initPSD_uni,c_uni,lambda_uni] = plotInitLeastErrorPSD(x_f,f_uni,xi_uni,M_uni,mMax_uni);
% figure(1)
% ylabel('$v_{\mathrm{L}}^{\mathrm{uni}}$ [mm$^2$]','FontSize',18,'Interpreter','Latex')
% xlabel("$L$ [mm]")
% stem(xi_uni,w_uni./(xi_uni), "filled", "DisplayName", "$\widetilde{v}_{\mathrm{L}}(t=0)/L_{\alpha}$", "LineWidth",2)


[R_bi,E_bi,M_bi_,initx_bi,initPSD_bi,c_bi,lambda_bi] = plotInitLeastErrorPSD(x_f,f_bi,xi_bi,M_bi,mMax_bi);
figure(1)
ylabel('$v_{\mathrm{L}}^{\mathrm{bi}}$ [mm$^2$]','FontSize',18,'Interpreter','Latex')
xlabel("$L$ [mm]")
stem(xi_bi,w_bi./(xi_bi), "filled", "DisplayName", "$\widetilde{v}_{\mathrm{L}}(t=0)/L_{\alpha}$", "LineWidth",2)

E_uni(E_uni==1) = [];
E_L_uni = sum(E_uni)/length(E_uni);
% E_bi(E_bi==1) = [];
% E_L_bi = sum(E_bi)/length(E_bi);