%% Adaptive Wheeler algorithm
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Implementation from:
% Computational Models for Polydisperse Particulate and Multiphase Systems,
% Marchisio and Fox (2013)
%
% compute quadrature approximation with the Wheeler algorithm (Sack &
% Donovan (1971)),
%
% INPUT:  N_delta       number of nodes of the quadrature approximation 
%         m             moments from 0 to 2N_delta-1 [m(1), ..., m(2N)] 
% 
% OUTPUT: xi_alpha      nodes of the N_delta disperse size distribution
%         w_alpha       weights of the N_delta disperse size distribution

function [xi_alpha,w_alpha]=Wheeler(m,N_delta) 

% Calculate intermediate quantities
for l=1:2*N_delta-2
    sigma(1,l+1) = 0.0;
end
%
for l=0:2*N_delta-1
    sigma(2,l+1) = m(l+1);
end
% 
% compute coefficients for Jacobi matrix 
%
a(1) = m(2)/m(1);
b(1) = 0.0;
%
for k=1:N_delta-1 
    for l=k:2*N_delta-k-1 
        sigma(k+2,l+1) = sigma(k+1,l+2)-a(k)*sigma(k+1,l+1)-b(k)*sigma(k,l+1); 
        a(k+1) = -sigma(k+1,k+1)/sigma(k+1,k)+sigma(k+2,k+2)/sigma(k+2,k+1);
        b(k+1) = sigma(k+2,k+1)/sigma(k+1,k);
    end
end
% 
% compute Jacobi matrix 
%
for i=1:N_delta 
    jacobi(i,i)=a(i);
end
%
for i=1:N_delta-1
    jacobi(i,i+1) = -(abs(b(i+1)))^0.5;
    jacobi(i+1,i) = -(abs(b(i+1)))^0.5;
end
% 
%compute eigenvalues and eigenvectors 
%
[evec,eval]=eig(jacobi); 
% 
% return weights 
%
for i=1:N_delta 
    w_(i)=evec(1,i)^2*m(1);
end
% 
% return abscissas 
%
for i=1:N_delta 
    xi_(i)=eval(i,i); 
end

w_alpha = w_';
xi_alpha = xi_';
end