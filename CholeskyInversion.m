%% Cholesky inversion of a Hessian matrix
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Checks if a Hessian matrix is positive definite, but does not determine other 
% definiteness (Indefiniteness, semidifiniteness). If the check is 
% successfull it will compute the inverse of the Hessian, else it will return
% zero.
% 
% INPUT:   Hessian       The Hessian matrix of a function

% OUTPUT:  H_inv         Inverse of the Hessian Matrix if positive
%                        definitene, else H_inv = 0
%          PosDef        integer flag that is zero if the Hessian is positive
%                        definite and if not positive definite it is a
%                        random positive integer

function [H_inv, PosDef] = CholeskyInversion(Hessian)
% Try Cholesky decomposition
[CholeskyMatrix,PosDef] = chol(Hessian);
% If positive definite chol(Hessian) will respond with zero and an inverse 
% can be computed
if ~PosDef
    H_inv = (CholeskyMatrix^(-1)) * CholeskyMatrix^(-1)';    
else
    H_inv = 0;
    return
end


