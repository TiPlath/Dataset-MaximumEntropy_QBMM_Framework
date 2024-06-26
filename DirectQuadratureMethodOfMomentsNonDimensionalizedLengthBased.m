%% Non-dimensionalised Direct Quadrature Method of Moments length-based
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Computes one time step of the direct quadrature method of moments 
% (According to Marchisio 2003) for aggregation, breakage and growth
% with a length-based kernel

% INPUT: xi_alpha       a vector containing the nodes of the discrete
%                       quadrature distribution
%        w_alpha        a vector containing the weights of the discrete
%                       quadrature distribution
%        dt             a double containing the discrete timestep
%        G              a vector containing the growth rate for each node
%        a              a vector containing the aggregation rate for each
%                       node
%        beta           a vector containing the breakage rate for each node
%        b_alpha        an anonymous function which computes the
%                       fragmentation distribution function for each node
%                       and moment
%        N_f            a double containing the number of fragments
%                       obtained from the fragmentation distribution
%                       function (zeroth moment)

% OUTPUT: xi            a vector containing the discrete quadrature distribution
%                       nodes of the next timestep
%         w             a vector containing the discrete quadrature distribution
%                       nodes of the next timestep

function [xi,w] = DirectQuadratureMethodOfMomentsNonDimensionalizedLengthBased(xi,w,dt,g,a,beta,b_alpha,N_f)
    n = length(w);
    %% source term computation and timestepping
    % Select a suitable kernel
    [A1,A2,S_L] = GrowthHydrodynamicAggregationPowerLawBreakageLengthBased(xi,w,g,a,beta,b_alpha,N_f);
    % solve equation for a, b
    A_ = [A1 A2];
    alp = A_\S_L;
    a_=alp(1:n);
    b=alp(n+1:end);
    % initialize the weighted abscissae
    s=w.*xi;
    % evolve weights, abscissae
    w=w+a_*dt;
    s=s+b*dt;
    xi = s./w;
end