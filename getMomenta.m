%% Get moments of a distribution split into dirac-delta distributions
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% takes a set of weights w_alpha and nodes xi_alpha and computes momenta m_k 
% such that the distribution 
%    f(xi) = sum_alpha^N_delta w_alpha delta(xi-xi_alpha) 
% has the requested moments, i.e. 
%    m_k = sum_alpha^N_delta w_alpha*xi_alpha^k
% 
% INPUT: xi_alpha      nodes of the discrete quadrature size distribution
%        w_alpha       weights of the discrete quadrature size distribution
% 
% OUTPUT: M            A vector with desired Moments of the current distribution
%                      with the range from 2...50 moments (1x50 at max)

function M = getMomenta(xi_alpha,w_alpha)
M = zeros(2*length(xi_alpha),1);
for k=0:2*length(w_alpha)-1
    M(k+1) = dot(w_alpha,xi_alpha.^k);
end
end
