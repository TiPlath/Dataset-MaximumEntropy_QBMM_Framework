%% Validate cumulative density functions
% Author = Plath, Timo
% E-mail: t.plath@utwente.nl
% Version = 1.0
%
% Checks if a cumulative density function starts with a zero
%
% INPUT:  CDF               matrix of a cumulative density function to validate
% 
% OUTPUT: validatedCDF      matrix of a validated cumulative density function

function [validatedCDF] = validateCDF(CDF)

% insert zero at the beginning
if CDF(1,2) ~= 0
    rowToInsert = [CDF(1,2),0];
    validatedCDF = [rowToInsert;CDF];
end
validatedCDF = CDF;
end