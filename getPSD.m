%% get Maximum Entropy particle size distribution
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

%Implementation from:
%Comparisons of methods for reconstructing particle size distribution from its moments
%by Shaohua Wu, Wenming Yang (2019)

% We use Lagrangian multipliers to find a distribution consisting of a series
% of exponential polynomials. We try to find the multipliers sucht that
% the exponential function has the same moments as our distribution to reconstruct

% INPUT: domain                  a vector containing the reconstruction domain
%        M                       a vector containing the Moments of the discrete
%                                distribution
%        mMax                    integer for the maximum number of moments

% OUTPUT: PSD                    Reconstructed PSD N(x) as a matrix with 
%                                different entries of calculations with a 
%                                different number of moments in every loop,
%                                in range of domain = xmin...xmax
%          PSDErrorComputation   a matrix of a PSD with a resolution
%                                matching the initial domain resolution,
%                                mainly used for error computation
%          k                     integer for the highest number of moments
%                                used in the reconstruction
%          E                     a vector of moment errors for each moment
%          lambda1               a vector of lambda values for the last
%                                reconstruction before exiting this
%                                function


function [PSD,PSDErrorComputation,k,E,lambda1] = getPSD(domain,M,mMax)
% number of points for the continuous reconstruction function
resolution = 300;
% Initialise PSD output matrix
PSD = zeros(mMax,resolution);
% Initialise PSD error computation matrix. This avoids having to 
% interpolate for the error computation.
PSDErrorComputation = zeros(mMax,length(domain));
% Initialize error output vector
E = ones(1,mMax);

for k = 2:mMax
    % zeroth moment
    M0 = M(1,1);
    % Determine accuracy
    tolerance = 10^(-4);  
    % change of variables to Shift interval from [xmin...xmax] to [0...1] since
    % Gaussian quadrature points are computed for normalised densities at
    % interval [0...1]
    Mur = zeros(1,k);
    for r = 0:k-1
        Mur(1,r+1) = (M(1,r+1))/((max(domain))^r*M0);
    end
    % iteratively fit all Lagrange multipliers to the corresponding moments
    for Z = 1:k-1
        % Initialize Hessian Matrix
        H = zeros(Z+1,Z+1);
        % Initialize vector for approximated moments
        I = zeros(Z+1,1);
        % Initialize vector for initial Lambda
        lambda0 = zeros(Z+1,1);
        % Initialize vector for corrected Lambda
        lambda1 = ones(Z+1,1);
        % Start Newton's iteration loop
        while norm(lambda1(1:Z,1) - lambda0(1:Z,1)) > tolerance
            % Set Parameter to iterate
            lambda0 = lambda1;
            % calculate approximated moments via Gauss quadrature (Jacobian)
            for n = 0:Z
                I = Gauss(I,Z,lambda0,n);
            end
            % calculate Hessian of lambda function via Gauss quadrature
            for m = 0:Z
                for n = 0:Z
                    H = Gauss(H,Z,lambda0,n,m);
                end
            end
            
            % Check if H is positive definite and invert via Cholesky decomposition
            [H_inv,PosDef] = CholeskyInversion(H);
            % If not positive definite Cholesky decompoisition returns
            % integer PosDef > 0
            if PosDef > 0
                % sprintf("Hessian matrix is not positive definite for %d Moments. Breaking Newton's Loop",k)
                PSD(k:mMax,:) = 0;
                k = k-1;
                return
            end
            % Update Lagrange multipliers with a damped newton solver to avoid
            % overshooting
            lambda1(1:Z+1,1) = lambda0(1:Z+1,1) - (H_inv*(Mur(1,1:Z+1)'-I));
        end
        E(1,1) = norm((Mur(1:length(I)) - I(1:end)')./Mur(1:length(I)));
    end
    % Relative moment error (Determines goodness of approach)
    E(1,k) = norm((Mur(1:length(I)) - I(1:end)')./Mur(1:length(I)));

    % increase domain resolution by linear interpolation
    domain_interp = interp1(linspace(0, 1, size(domain, 2)), domain, linspace(0, 1, resolution), 'linear', 0);
    % define anonymous function for the continuously reconstructed PSD
    PSDFunction = @(x) (exp(-(sum(lambda0 .* (x / max(x)).^((0:length(lambda0)-1)'))))) * M0 / max(x);
    % assign reconstructed PSD with a certain resolution of bins to the PSD output vector
    PSD(k,:) = PSDFunction(linspace(min(domain_interp),max(domain_interp),length(domain_interp)));
    PSDErrorComputation(k,:) = PSDFunction(domain);
end
end


