%% Non-dimensionalised Direct Quadrature Method of Moments
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Computes one time step of the direct quadrature method of moments 
% (According to Marchisio 2003) for aggregation, breakage and growth
% with a volume-based kernel

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
%                       

% OUTPUT: xi_alpha      a vector containing the discrete quadrature distribution
%                       nodes of the next timestep
%         w_alpha       a vector containing the discrete quadrature distribution
%                       nodes of the next timestep

function [xi_alpha,w_alpha] = DirectQuadratureMethodOfMomentsNonDimensionalized(xi_alpha,w_alpha,dt,g,a,beta,b_alpha,N_f)
    n = length(w_alpha);
    %% source term computation and timestepping
    % Select a suitable kernel
    [A1,A2,S_V] = GrowthAggregationBreakageNonDimensionalized(xi_alpha,w_alpha,g,a,beta,b_alpha,N_f);
    % solve equation for a, b
    A_ = [A1 A2];
    alp = A_\S_V;
    a_=alp(1:n);
    b=alp(n+1:end);
    % initialize the weighted abscissae
    s=w_alpha.*xi_alpha;
    % evolve weights, abscissae
    w_alpha=w_alpha+a_*dt;
    s=s+b*dt;
    xi_alpha = s./w_alpha;
end