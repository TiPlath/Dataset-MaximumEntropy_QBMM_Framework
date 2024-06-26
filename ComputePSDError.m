%% Compute the reconstruction Jenson-Shannon divergence for a fixed c-value

% We perform a maximum entropy reconstruction for a given c-value which
% determines the maximum particle size. For the reconstructed distributions
% with 2...mMax moments, the Jenson-Shannon divergence is computed against the
% initial distribution to obtain the reconstruction with the least error.
% This function mainly serves the bisection method as input to find the
% optimal c-value
% 

% INPUT: f           a vector containing the initial density distribution values
%        x_f         a vector containing the initial density distribution classes
%        c           a double containing the multiplier for the relative
%                    distance between the largest node and the maximum particle 
%                    size of the reconstruction domain
%                    reconstruction domain
%        xi          a vector containing the nodes of the discrete quadrature
%                    distribution
%        M           a vector containing the moments of the initial
%                    distribution
%        mMax        an integer containing the maximum number of moments

% OUTPUT: PSD        a vector containing the reconstructed PSD with the least 
%                    Jenson-Shannon divergence compared to the initial distribution
%         x_interp   a vector containing the reconstruction domain which
%                    is linearly interpolated to match the resolution of
%                    the PSD
%         RMin       a double containing the minimum Jenson-Shannon divergence
%                    for the current c-value
%         RPos       a scalar containing the position (= number of moments) 
%                    of the minimum Jenson-Shannon divergence

function [PSD,x_interp,E,RMin,RPos,lambda] = ComputePSDError(x_f,f,c,xi,M,mMax)
    % calculate maximum particle size
    xMax = c*max(xi);
    xMin = 0;
    % Global domain with respect to the size of particles
    x = linspace(xMin, xMax, length(f));
    % Compute reconstructed distribution for different numbers of moments
    [PSD,PSDErrorComputation,k,E,lambda] = getPSD(x, M, mMax);
    % initialise error vector
    R = zeros(1,mMax);    
    % intialise interpolated PSD to compute Jenson-Shannon divergence
    PSDInterpolated = zeros(mMax,length(x));
    % Estimation error between real and approximated PSD by Jenson-Shannon
    % divergence
    JSD = zeros(1,mMax);
    for i = 1:mMax
        % Interpolate Grid x to x_f for error computation.
        PSDInterpolated(i,:) = interp1(x,PSDErrorComputation(i,:),x_f, 'linear', 0);
        JSD(i) = jensenShannonDivergence(f(1,:),PSDInterpolated(i,:));
    end
    % interpolated domain to the resolution of PSD. This is mainly done to
    % make the plots look smoother
    x_interp =  interp1(linspace(0, 1, size(x, 2)), x, linspace(0, 1, length(PSD)), 'linear', 0);
    % find minimum of R
    [RMin,RPos] = min(JSD);
end