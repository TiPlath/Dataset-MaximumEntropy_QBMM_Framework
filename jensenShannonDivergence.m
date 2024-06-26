%% Computes the Jensen-Shannon Divergence between two probability distributions.
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0

% Computes the Jensen-Shannon Divergence between two probability distributions
% P and Q. P and Q should be vectors representing probability distributions.

% The Jensen-Shannon Divergence is a measure of similarity between two distributions.
% It is defined as half the sum of the Kullback-Leibler divergences between P
% and the average distribution M of P and Q, and between Q and M.

% If either P or Q is a zero vector (indicating an empty distribution),
% the Jensen-Shannon Divergence is set to 1, representing maximum dissimilarity.

% INPUT:  P,Q       vector of density distribution function values
% 
% OUTPUT: JSD       Jensen-Shannon divergence measuring the similarity of P
%                   and Q in range 0,1.


function JSD = jensenShannonDivergence(P, Q)

    % Ensure that distributions are normalised
    if (sum(Q) == 0 || sum(P) == 0)
        JSD = 1; % Set JSD to 1 for maximum dissimilarity if one distribution is empty
        return
    end
    P = P / sum(P);
    Q = Q / sum(Q);
    
    % Calculate the average distribution M
    M = 0.5 * (P + Q);
    
    % Calculate the KL divergence of P from M
    D_KL_P_M = sum(P .* log2(P ./ M), 'omitnan');
    
    % Calculate the KL divergence of Q from M
    D_KL_Q_M = sum(Q .* log2(Q ./ M), 'omitnan');
    
    % Calculate the Jensen-Shannon Divergence
    JSD = 0.5 * (D_KL_P_M + D_KL_Q_M);
end

