%% Compute kernels for constant growth, aggregation and breakage
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

function [A1,A2,S_V] = GrowthAggregationBreakageNonDimensionalized(xi,w,G,a,beta,b_alpha,N_f)
    % number of weights
    n = length(w);
    %% Assemble source term
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
        A1(k+1,:)= (1-k)*xi.^k;
        A2(k+1,:)= k*xi.^(k-1);
        for i = 1:n
            for r = 1:n
                % volume-based aggregation
                A3(k+1,i) = A3(k+1,i) + (((xi(i) + xi(r))^k) * w(i) * w(r));
                A4(k+1,i) = A4(k+1,i) + ((xi(i)^k) * w(i) * w(r));
            end
        end
        % volume-based breakage
        A5(k+1,:) = xi.^k;
        B1(k+1,:)= b_alpha(xi,k,N_f);   
    end 
    %% compute source terms
    % Aggregation
    Ag = 0.5*(A3*a) - A4*a;
    % Breakage 
    Br = B1*diag(w)*beta - A5*diag(w)*beta;
    % Growth
    Gr = A2*diag(w)*G;
    % breakage, aggregation and growth
    S_V = Ag + Br + Gr;
end