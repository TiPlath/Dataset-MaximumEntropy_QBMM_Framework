%% Gauss Quadrature with 24 points.
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% A Gauss-Legendre quadrature is employed to compute either the gradient of
% our Lagrangian functional or Hessian matrix depending on the number of
% input arguemnts.

% INPUT: f              a preallocated matrix with entries for the gradient
%                       or Hessian matrix to compute
%        Z              an integer which determines the number of
%                       Lagrangian multipliers
%        lambda0        a vector of Lagrangian multipliers
%        n              an integer which determines the power for the
%                       matrix entry to compute 
%        m              an integer which determines the power for the
%                       matrix entry to compute 

% OUTPUT: f             a matrix with a computed value for a certain entry
%                       determined by n for the gradient and nxm for the
%                       Hessian matrix


function f = Gauss(f,Z,lambda0,n,m)

% Gauss nodes
xi = [-0.064056893, 0.064056893, -0.191118867, 0.191118867, -0.31504268, ...
      0.31504268, -0.433793508, 0.433793508, -0.545421471, 0.545421471, ...
      -0.648093652, 0.648093652, -0.740124192, 0.740124192, -0.820001986, ...
      0.820001986, -0.886415527, 0.886415527, -0.938274552, 0.938274552, ...
      -0.974728556, 0.974728556, -0.99518722, 0.99518722];
  
% Weights
wi = [0.127938195, 0.127938195, 0.125837456, 0.125837456, 0.121670473, ...
      0.121670473, 0.115505668, 0.115505668, 0.1074427, 0.1074427, 0.097618652, ...
      0.097618652, 0.086190162, 0.086190162, 0.073346481, 0.073346481, ...
      0.059298585, 0.059298585, 0.044277439, 0.044277439, 0.028531389, ...
      0.028531389, 0.01234123, 0.01234123];

% Computing the approximated moments by Gaussian quadrature
if nargin == 4
    ISum = 0;
    for i = 1:length(xi)
        expSum = 0;
        for j = 0:Z
            expSum = expSum + lambda0(j+1,1)*(0.5*xi(i) + 0.5)^j;
        end        
        ISum = ISum + wi(i)*(0.5*xi(i)+0.5)^n * exp(-expSum); 
    end
    f(n+1) = 0.5 .* ISum;
end

% Computing the Hessian matrix by Gaussian quadrature
if nargin == 5
    HSum = 0;
    for i = 1:length(xi)
        expSum = 0;
        for j = 0:Z
            expSum = expSum + lambda0(j+1,1)*(0.5*xi(i)+0.5)^j;
        end
        HSum = HSum + wi(i)*(0.5*xi(i)+0.5)^(n+m) * exp(-expSum); 
    end
    f(m+1,n+1) = 0.5 .* HSum;
end
    
    