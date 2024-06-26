%% Convert a PDF to a CDF
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Convert a cumulative density function (CDF) to a (probability) density 
% function (PDF) using a right Riemann sum.
% 
% INPUT: CDF           a matrix containing the cumulative density
%                      distribution function values and classes
% 
% OUTPUT: PDF          a matrix containing the density distribution
%                      function values and classes

function PDF = convertCDFtoPDF(CDF)
    PDF = zeros(length(CDF),2);
    PDF(:,1) = CDF(:,1);
    PDF(:,2) = [diff(CDF(:,2))./diff(CDF(:,1)); 0];