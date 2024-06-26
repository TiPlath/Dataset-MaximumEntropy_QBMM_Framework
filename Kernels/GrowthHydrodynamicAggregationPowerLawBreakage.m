%% Compute kernels for constant growth, hydrodynamic aggregation and power law breakage
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Computes the source term for constant aggregation, constant breakage
% and constant growth with a volume-based kernel

% INPUT: xi_alpha       a vector containing the nodes of the discrete
%                       quadrature distribution
%        w_alpha        a vector containing the weights of the discrete
%                       quadrature distribution
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

% OUTPUT: A1            a matrix containing the prefactors of the 
%                       time-dependent weight term for each node and moment
%         A2            a matrix containing the prefactors of the
%                       time-dependent weighted node term for each node and
%                       moment
%         S_V           a vector containing the volume-based source term 
%                       contributions to each of the moments, i.e. dM_k/dt = S_V

function [A1,A2,S_V] = GrowthHydrodynamicAggregationPowerLawBreakage(xi_alpha,w_alpha,dotV,a,beta,b_alpha,N_f)
    % number of weights
    n = length(w_alpha);
    %% Assemble source term
    I = ones(n,1);
    % assemble matrices A1, A2 (left hand side of equation)
    A1=zeros(2*n,n);
    A2=zeros(2*n,n);
    % specific for aggregation (right hand side of equation)
    A3=zeros(2*n,n);
    A4=zeros(2*n,n);
    % specific for breakage (right hand side of equation)
    A5=zeros(2*n,n);
    % fragmentation distribution function
    B1=zeros(2*n,n);
    for k=0:2*n-1
        % volume-based growth
        A1(k+1,:)= (1-k)*xi_alpha.^k;
        A2(k+1,:)= k*xi_alpha.^(k-1);
        % volume-based hydrodynamic aggregation
        for i = 1:n
            for r = 1:n
                A3(k+1,i) = A3(k+1,i) + ((xi_alpha(i) + xi_alpha(r))^k)*(a(i)*(xi_alpha(i)^(1/3)+xi_alpha(r)^(1/3))^(3)) * w_alpha(i) * w_alpha(r);
                A4(k+1,i) = A4(k+1,i) + (xi_alpha(i)^k) *(a(i)*(xi_alpha(i)^(1/3)+xi_alpha(r)^(1/3))^(3)) * w_alpha(i) * w_alpha(r);
            end
        end
        % volume-based power law breakage
        A5(k+1,:) = beta.*xi_alpha.^(1/3).*xi_alpha.^k;
        B1(k+1,:) = beta.*xi_alpha.^(1/3).*b_alpha(xi_alpha,k,N_f);
    end
    %% compute source terms
    % Aggregation
    Ag = 0.5*(A3)*I - A4*I;
    % Breakage 
    Br = B1*diag(w_alpha)*I - A5*diag(w_alpha)*I;
    % Growth (constant length growth)
    Gr =  A2*diag(w_alpha)*(dotV.*xi_alpha.^(1/3));
    % merge breakage, aggregation and growth contribution to source term
    S_V = Ag + Br + Gr;
end