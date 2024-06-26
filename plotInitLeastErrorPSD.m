%% Initial reconstruction bisection approach to determine an optimal c-value
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Reconstructs the PSD by a maximum entropy approach and gets the PSD with
% the least error accounting for maximum particle size and highest possible
% number of moments for the initial time step. We employ a bisectional
% approach to determine the optimal c-value which has the least error.

% INPUT:    x_f        a vector containing the domain of the measured PSD
%           f          a vector containing the density values of the measured PSD
%           xi         a vector which contains the nodes from Wheeler of the 
%                      real PSD
%           M          a Matrix which contains the moments of the real PSD
%           mMax       The maximum number of moments which where computed

% OUTPUT:   R          Jenson-Shannon divergence determining the equality
%                      of two distributions
%           E          a vector of moment errors of the maximum entropy reconstruction
%           M          a vector with the least error moments
%           x          a vector with the domain which gives the least error for the
%                      initial time step
%           PSD        a vector which gives the reconstructed discrete particle size
%                      distribution with the least error against the real PSD
%           c          a value which holds the relative distance of the
%                      largest node to the maximum particle size of the 
%                      reconstruction domain
%           lambda     a vector of obtained lambdas from the reconstruction
%                      function
function [R,E,M,x,PSD,c,lambda] = plotInitLeastErrorPSD(x_f,f,xi,M,mMax)

%%%%%%%%%%%%%%%%% Set initial parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize factor c which will help to determine the maximum particle size
cLow = 1;
cHigh = 3;
c = (cLow + cHigh)/2;
% Number of iteratisons (very sensible to this)
it = 5;
% Accuracy of the iteration as alternative stopping criterion
tolerance = 10^(-20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start while loop to iteratively determine c to calculate maximum particle size
for iterations = 1:it
    [PSD1,x1,E1,RMin1,RPos1,lambda1] = ComputePSDError(x_f',f',c,xi,M,mMax);
    [PSD2,x2,E2,RMin2,RPos2,lambda2] = ComputePSDError(x_f',f',cHigh,xi,M,mMax);
    if RMin1 < tolerance || iterations == it
        break
    elseif RMin1 > RMin2
        cLow = c;
    else
        cHigh = c;
    end
    if iterations == 5
        break
    end
    c = (cLow + cHigh)/2;
end
% Allocation of variables with least error
if RMin1 < RMin2
    PSD = PSD1(RPos1,:)';
    M = M(1:RPos1);
    x = x1';
    E = E1;
    R = RMin1;
    c = c;
    lambda = lambda1;
else
    PSD = PSD2(RPos2,:)';
    M = M(1:RPos2);
    x = x2';
    E = E2;
    R = RMin2;
    c = cHigh;
    lambda = lambda2;
end
%% Plot reconstructed PSD against Analytical solution
figure(1)
axes1 = axes('Parent', figure(1));
plot(x_f,f, 'LineWidth', 2);
xlim([0 max(x_f)])
legend('Analytic', 'location', 'best')
hold on
plot(x(1:end),PSD(1:end), '--', 'LineWidth', 1.5, 'Marker', 'x', 'MarkerIndices', 1:4:length(x));
xlabel('$L$','FontSize',18,'Interpreter','Latex')
ylabel('$v_{\mathrm{L}}$','FontSize',18,'Interpreter','Latex')
legend('$v_{\mathrm{L}}(t = 0)$', '$\bar{v}_{\mathrm{L}}(t = 0)$', 'location', 'best','Interpreter','Latex')
set(gcf, 'Color', 'w');
set(axes1,'FontSize',16)