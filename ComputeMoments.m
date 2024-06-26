%% Compute moments of a density distribution
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Computes the moments of a given density distribution and returns them in
% a vector M. The computed moments are raw moments according to the formula
% M_k = \int_x (x^r * f(x))
% 
% INPUT:   f             a vector containing the density distribution values
%          x             a vector containing the density distribution classes
%          mMax          an integer containing the maximum number of moments 
% 
% OUTPUT:  M             a vector with desired moments of the current density
%                        distribution

function M = ComputeMoments(x, f, mMax)


% Initialize moment matrix
M = zeros(1,mMax);
for r = 0:mMax-1
    M(r+1) = sum(diff(x).*f(1:end-1).*(x(1:end-1).^r));   
end
end