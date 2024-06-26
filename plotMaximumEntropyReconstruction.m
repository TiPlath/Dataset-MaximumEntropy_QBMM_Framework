%% Plot the maximum entropy reconstruction
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Gets the reconstructed PSD from getPSD() and plots the reconstructed distribution 
% with the least moment error and makes sure that the number of moments is
% uneven to better represent skewed distributions.
% 
% INPUT:  xi            node vector from Wheeler
%         w             weight vector from Wheeler
%         mMax          maximum number of moments
%         dxMax         relative distance of the reconstruction domain to
%                       the largest node L_{N_{\delta}}
%         f             initial distribution function, used to obtain the
%                       length of the discretised reconstruction domain
% 
% OUTPUT: PSD           reconstructed particle size distribution
%          x_interp     reconstruction domain interpolated to the
%                       continuous reconstruction resolution
%          k            number of moments taken for reconstruction
%          lambda       lambda values of the reconstruction function

function [PSD,x_interp,k,lambda] = plotMaximumEntropyReconstruction(xi,w,mMax,xMin,dxMax,f)
    
    % get momenta of PSD
    m(1,:) = getMomenta(xi,w);
    
    % Set domainaperpae
    x = linspace(xMin,max(xi)+dxMax,length(f));
    %% get PSD in a maximum entropy reconstruction
    [PSD,~,k,E,lambda] = getPSD(x, m, mMax);
    % choose lowest moment error PSD
    [~,EPos] = min(E(2:end));
    % correct EPos by +1 due to ignoring the first entry
    EPos = EPos + 1;
        if mod(EPos,2)==0 && EPos ~= 2
            EPos = EPos-1;
            k = k-1;
        end
    PSD = PSD(EPos,:);
    % interpolate domain to the resolution of PSD. This is mainly done to
    % make the plots look smoother
    x_interp =  interp1(linspace(0, 1, size(x, 2)), x, linspace(0, 1, length(PSD)), 'linear');
    %% Plot reconstructed PSD with least moment error and uneven moment number
    figure(1)
    plot(x_interp,PSD, '--', 'LineWidth', 1)
    drawnow
end