%% Convert a PDF to a CDF
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Convert a (probability) density function (PDF) to a cumulative density
% function (CDF) using a right Riemann sum.
% 
% INPUT: PDF           a matrix containing the density distribution
%                      function values and classes
% 
% OUTPUT: CDF          a matrix containing the cumulative density
%                      distribution function values and classes

function CDF = convertPDFtoCDF(PDF)
    CDF = zeros(length(PDF),2);
    CDF(:,1) = PDF(:,1);
    for i = 2:length(PDF)
        CDF(i,2) = CDF(i-1,2) + (PDF(i,1) - PDF(i-1,1)) * PDF(i-1,2);
    end