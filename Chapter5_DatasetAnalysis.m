%% Application to twin-screw wet granulation data: Dataset analysis (1) (section 5.2)

% We perform a dataset analysis of the twin-screw wet granulation dataset 
% (Plath et al. 2021) to obtain functions that qualitatively match the experimental
% results of the processed lactose, polyvinylpyrrolidon and microcrystalline
% cellulose mixture.
% 1. We iteratively determine, the aggregation rate prefactor 
%    a_{t,L} such that steady state is reached in the mean residence time.
% 2. Keeping the aggregation rate prefactor constant, the 
%    breakage rate prefactor beta_{t,L}' is set to match the 
%    Sauter mean diameter $d_{32}$.
% 3. Since the breakage rate also influences the time it takes to reach 
%    steady state, one has to readjust the aggregation rate prefactor. Both 
%    prefactors are set to match MRT and d_32 with a relative error of 0.1%.

clear all
close all
clc
% turn off warning for nearlySingularMatrices which occur naturally in DQMOM
warning_state = warning('off', 'MATLAB:nearlySingularMatrix');
% add the functions in the kernel path to the script
addpath(genpath('Kernels'));
%% define variables (pre-processing)
% try to get the aggregation and breakage that fits the Sauter mean
% diameter and the mean residence time of the respective distributions.
% number of particles
N_P = 8.78888e+07;
% dimensionless aggregation rate
aPrime = 1;
% parameters
% specific feed load values
SFL = [0.12,0.04,0.12,0.04,0.08];
% liquid to solid ratio values
LS = [0.35, 0.35, 0.1, 0.1, 0.225];
% assemble parameterPairs vector
parameterPairs = strings(1, length(SFL));
for i = 1:length(SFL)
    parameterPairs(i) = sprintf('$\\mathit{SFL}=%.2f$, $\\mathit{L/S}=%.2f$', SFL(i), LS(i));
end
% experimental points
numberOfExperiments = 5;
steadyStateTimes = zeros(numberOfExperiments,1);
Omega_t = zeros(numberOfExperiments,1);
realSimulationTime = zeros(numberOfExperiments,1);
% Sauter mean diameter
d32 = zeros(numberOfExperiments,1);
%% variables to fit
a_t = [2.97e4, 4.32e4, 6.41e4, 1.001e5, 5.1e4];
beta_tPrime = [0.0365,0.0508,0.3251,0.1995,0.1275];
MRT = [34.5,24.6,23.2,13.3,23.9];
%% Run different L/S and SFL ratios
% Aggregation rate is fit to match the steady state at MRT; Breakage rate is
% fit to match the Sauter mean diameter
experiment = {"MCCLactosePVP_HighSFL_HighLS_N16","MCCLactosePVP_LowSFL_HighLS_N14", ...
    "MCCLactosePVP_HighSFL_LowLS_N15","MCCLactosePVP_LowSFL_LowLS_N13","MCCLactosePVP_MediumSFL_MediumLS_N26"};
% define data storage
data = cell(1,numberOfExperiments);
% Loop over each experiment with its respective MRT
for i = 1:numberOfExperiments
    data{i} = struct('filename', experiment{i}, 'aggregationRate', a_t(i), ...
        'breakageRatePrime', beta_tPrime(i), "MRT", MRT(i), "steadyStateTimes", steadyStateTimes(i), ...
        "realSimulationTime", 0, "sauterDiameter", 0, "Omega_t", 0, "Omega_V", 0);
    fprintf("<strong>Simulating experiment %s and MRT = %.1fs</strong>\n", data{i}.filename, data{i}.MRT)
    [data{i}.steadyStateTimes,data{i}.realSimulationTime, data{i}.sauterDiameter, data{i}.Omega_t, data{i}.Omega_V]=Chapter5_MRTd32Match(data{i}.filename, data{i}.aggregationRate, data{i}.breakageRatePrime, data{i}.MRT);
    % save values to their respective vectors
    steadyStateTimes(i) = data{i}.steadyStateTimes;
    Omega_t(i) = data{i}.Omega_t;
    Omega_V(i) = data{i}.Omega_V;
    realSimulationTime(i) = data{i}.realSimulationTime;
    d32(i) = data{i}.sauterDiameter;
    breakageRate = beta_tPrime(1:numberOfExperiments)./Omega_t';
end
% non-dimensionalise d32
d32Prime = d32.*(Omega_V'.^(-1/3));
% aggregation rates in units per second
aggregationRateSeconds = a_t(1:numberOfExperiments).*N_P.*Omega_V;

%% Fit d32 against process parameter pairs SFL and LS.
rho_P = 1540;
bilinear = @(params, xy) (params(1).*xy(:,2).^2) .* (xy(:,1)./(1000*rho_P)).^(1/3);
px = SFL';
py = LS';
pz = d32;

guess = [1];
params = lsqcurvefit(@(params, xy) bilinear(params, xy), guess, [SFL', LS'], d32);

% Data for predicted vs measured plot
predicted = bilinear(params, [SFL', LS']);

% Compute R2
residuals = d32 - predicted;
ss_res = sum(residuals.^2);
ss_tot = sum((d32 - mean(d32)).^2);
r_squared = 1 - (ss_res / ss_tot);

figure(4);
ax = axes('Parent', figure(4));
plot(predicted, d32, 'k*', 'MarkerSize', 6,  'LineWidth', 2);
title(['R^2 = ' num2str(r_squared, '%.2f')]);
xlabel('Predicted', 'FontSize', 14,'Interpreter','Latex');
ylabel('Measured', 'FontSize', 14,'Interpreter','Latex');
grid on;
ax.FontSize = 14;

x = linspace(min(xlim), max(xlim)+0.1*max(xlim));
axis tight
hold on;
plot(x, x, 'k:', 'LineWidth', 2);
hold off;
set(gcf, 'Color', 'w');
set(ax,'FontSize', 16)

%% Contour plot d32 fit
            
[xData, yData, zData] = prepareSurfaceData( SFL, LS, d32 );

% Set up fittype and options.
ft = fittype( '(a*y^2)*(x/1.540)^(0.33333)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 1;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Make contour plot.
figure(5);
plot(fitresult, [xData, yData], zData, 'Style', 'Contour' );
% Label axes
xlabel('SFL', 'FontSize', 14, 'Interpreter', 'Latex');
ylabel('L/S', 'FontSize', 14, 'Interpreter', 'Latex');
grid off
title('$d_{32}$ Contour Plot', 'Interpreter', 'Latex');

% Customize colorbar
c = colorbar;
c.Label.String = "$d_{32}$ [m]";
c.Label.Interpreter = "Latex";
c.FontSize = 14;
c.Label.FontSize = 18;

% Customize axes
set(gca, 'FontSize', 14);

%% Fit aggregation rate against steady state time
syms aggregationRate(sym_MRT) sym_MRT
agg = @(a_c,p_a,t0,t) a_c.*((t-t0)./t0).^p_a;
a_c = 0.1879;
p_a = -0.6736;
t0 = 7.5;
aggregationRate(sym_MRT) = agg(a_c,p_a,t0,sym_MRT);

figure(6)
axes6 = axes('Parent', figure(6));
fplot(aggregationRate(sym_MRT), [10 50], 'LineWidth', 2, "DisplayName", "$a_{s}(\mathit{MRT}) = a_c \left(\frac{\mathit{MRT}-t_{opt}}{t_{opt}}\right)^{p_a}$")
hold on
legend('location','northeast', 'FontSize', 14,'Interpreter','Latex')
for j=1:numberOfExperiments
    plot(realSimulationTime(j),aggregationRateSeconds(j),"x", "DisplayName", parameterPairs(j),'MarkerSize', 14, 'LineWidth',2)
end
ylabel('$a_{\mathrm{s}}$ [s\textsuperscript{-1}]','FontSize', 14, 'Interpreter','Latex')
xlabel('$\mathit{MRT}$ [s]','FontSize', 14, 'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes6,'FontSize', 16)

%% Fit breakage rate against Sauter mean diameter
syms breakageRate_tPrime(sym_d32) sym_d32_2
beta_alpha = @(beta_c,p_b,x) beta_c.*x.^p_b;
beta_c = 2.62558;
p_b= -0.91461;
breakageRate_tPrime(sym_d32_2) = beta_alpha(beta_c,p_b,sym_d32_2);

figure(7)
axes7 = axes('Parent', figure(7));
fplot(breakageRate_tPrime(sym_d32_2), [5 120], 'LineWidth', 2, "DisplayName", "$\beta^{\prime}_{t}(d_{32}^{\prime}) = \beta_c {d_{32}^{\prime}}^{p_b}$")
hold on
ylim([0 0.6])
legend('location','northeast', 'FontSize', 14,'Interpreter','Latex')
for j=1:numberOfExperiments
    plot(d32Prime(j),beta_tPrime(j),"x", "DisplayName", parameterPairs(j),'MarkerSize', 14, 'LineWidth',2)
end
ylabel('$\beta^{\prime}_{\mathrm{t}}$ [-]','FontSize', 14, 'Interpreter','Latex')
xlabel('$d^{\prime}_{32}$','FontSize', 14, 'Interpreter','Latex')

set(gcf, 'Color', 'w');
set(axes7,'FontSize', 16)

%% plot zeroth qualityfactor over d_32
% qualityfactor of first fit
q0 = sqrt((aggregationRateSeconds./double(subs(aggregationRate(sym_MRT), realSimulationTime'))).^2);

residual_q0 = sum((q0 - 1).^2);
% Display residual q0
fprintf('Residual of q0: %.4f\n', residual_q0);

syms fit1(sym_d32)
fit = @(a,b,x) a.*x.^b;
a = 0.3478;
b = -0.1488;
fit1(sym_d32) = fit(a,b,sym_d32);

figure(8)
axes8 = axes('Parent', figure(8));
fplot(fit1(sym_d32), [1e-4 5e-3], 'LineWidth', 2, "DisplayName", "$a_{s}(\mathit{MRT}) = a_c \left(\frac{\mathit{MRT}-t_{opt}}{t_{opt}}\right)^{p_a}$")
hold on

legend('location','northeast', 'FontSize', 14,'Interpreter','Latex')
for j=1:numberOfExperiments
    plot(d32(j),q0(j),"x", "DisplayName", parameterPairs(j),'MarkerSize', 14, 'LineWidth',2)
end
ylabel('$\frac{a_{\mathrm{s}}}{a_{\mathrm{s}}(\mathit{MRT})}$ [-]','FontSize', 14, 'Interpreter','Latex')
xlabel('$d_{32}$ [m]','FontSize', 14, 'Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes6,'FontSize', 16)

%% plot quality factors in a single plot

figure(9)
hold on

q1 = sqrt((q0./double(subs(fit1(sym_d32),d32'))).^2);

residual_q1 = sum((q1 - 1).^2);
% Display residual q1
fprintf('Residual of q1: %.4f\n', residual_q1);

% fplot(fit1(sym_MRT), [5 50])
plot(d32,q1, 'x', 'MarkerSize', 14, "DisplayName", "$q_1$")
plot(d32,q0, '+', 'MarkerSize', 14, "DisplayName", "$q_0$")
legend('Interpreter','Latex')

%% Aggregation rate contour plot
% Prepare data
[xData, yData, zData] = prepareSurfaceData(MRT, d32, aggregationRateSeconds);

% Define the custom power-law model
powerLawModel = @(coeffs, xy) coeffs(1) * ((xy(:,1)-7.5)/7.5).^coeffs(2) .* (xy(:,2)/Omega_V(1)^(1/3)).^coeffs(3);

% Initial parameter guess
initialParams = [0.186872604554379, 0.489764395788231, 0.445586200710899];

% Fit model to data using lsqcurvefit
fitresult = lsqcurvefit(powerLawModel, initialParams, [xData, yData], zData);

% Display the results
disp(fitresult);

% the fitted parameters
params = [fitresult(1),fitresult(2),fitresult(3)];

% Generate a grid of values for x and y
xValues = linspace(min(xData), max(xData), 50);
yValues = linspace(min(yData), max(yData), 50);
[X, Y] = meshgrid(xValues, yValues);

% Evaluate the model at each point in the grid
Z = powerLawModel(params, [X(:), Y(:)]);

% Reshape the result to match the grid
Z = reshape(Z, size(X));

% Create the contour plot
ax2 = figure(10);
contourf(X, Y, Z,"EdgeColor", "none", 'LevelStep',0.03);
hold on
plot(xData,yData, 'o', 'MarkerSize', 9, 'MarkerFaceColor',[1 1 1], 'MarkerEdgeColor',[0 0 0]);
% Label axes
xlabel('$\mathit{MRT}$ [s]', 'Interpreter', 'Latex');
ylabel('$d_{32}$ [m]', 'Interpreter', 'Latex');

c = colorbar;
c.Label.String = "$a_{\mathrm{s}}$ [s\textsuperscript{-1}]";
c.Label.Interpreter = "Latex";
c.FontSize = 14;
c.Label.FontSize = 18;

ax2.CurrentAxes.FontSize = 14;

%% Fit a_s against d_32 and MRT measured vs. predicted.

% Data for predicted vs measured plot
predicted = powerLawModel(params, [MRT', d32]);

% Compute R2
residuals = aggregationRateSeconds' - predicted;
ss_res = sum(residuals.^2);
ss_tot = sum((aggregationRateSeconds - mean(aggregationRateSeconds)).^2);
r_squared = 1 - (ss_res / ss_tot);

ax2 = figure(11);
ax2 = axes('Parent', figure(11));
plot(predicted, aggregationRateSeconds, 'k*', 'MarkerSize', 6,  'LineWidth', 2);
title(['R^2 = ' num2str(r_squared, '%.3f')]);
xlabel('Predicted', 'FontSize', 14,'Interpreter','Latex');
ylabel('Measured', 'FontSize', 14,'Interpreter','Latex');
grid on;
ax2.FontSize = 14;

x = linspace(min(xlim), max(xlim)+0.1*max(xlim));
axis tight
hold on;
plot(x, x, 'k:', 'LineWidth', 2);
hold off;

set(gcf, 'Color', 'w');
set(ax2,'FontSize', 16)