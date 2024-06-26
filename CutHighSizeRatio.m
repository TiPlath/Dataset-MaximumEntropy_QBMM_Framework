function cutPDF = CutHighSizeRatio(PDF,quantileMin,quantileMax,minPolydispersity)

if quantileMin == 0 && quantileMax == 1
    cutPDF = PDF;
    return;
end
% maximum probability
maxProb = PDF(end,2);
% normalize
PDF(:,2) = PDF(:,2)/PDF(end,2);
% Size ratio = rMax/rMin
SR = PDF(end,1)/PDF(1,1);
% fprintf('Initial size ratio is %f', SR)
radiusMin = getRadiusByQuantile(PDF,quantileMin);
radiusMax = getRadiusByQuantile(PDF,quantileMax);
% minimum polyidispersity at the base
radiusMinCut = min(radiusMin*(1+minPolydispersity),radiusMax);
% cut off min
while PDF(1,1) <= radiusMinCut
    PDF(1,:) = [];
end
PDF = [radiusMinCut,quantileMin;PDF];
if quantileMin ~= 0
    PDF = [radiusMin,0;PDF];
end
% cut off max
while PDF(end,1) >= radiusMax
    PDF(end,:) = [];
end
radiusMaxCut = max(radiusMax-(1-minPolydispersity)*(radiusMax-PDF(end,1)),radiusMin);
% PDF = [PDF;radiusMaxCut,quantileMax];
% linearly interpolate values to ensure non sharp cutoff
% dp1 = (PDF(end,2) - PDF(end-1,2));
% dp2 = 1 - PDF(end,2);
% dp = dp2/dp1;
% dp = round(dp);
% dx = (PDF(end,1)-PDF(end-1,1))/dp;

% for i = 1:dp
%     PDF = [PDF;PDF(end,1)+dx,PDF(end,2)+dp1];
% end
% delete columns above probability 1.0
while PDF(end,2) >= 1.0
    PDF(end,:) = [];
end
% delete end if radius above radiusMax
if PDF(end,1) >= radiusMax
    PDF(end,:) = [];
end
% PDF = [PDF;radiusMax,1];

SR = PDF(end,1)/PDF(1,1);
% fprintf('Cut size ratio is %f\n', SR)
%de-normalize to not loose information
PDF(:,2) = PDF(:,2)*maxProb;
cutPDF = PDF;
end




function radius = getRadiusByQuantile(PDF,quantile)
if quantile > 1 || quantile < 0
    error('quantile has to be between 0 and 1')
end
%Get readius index nearest to quantile
[~,nearestIndex] = min(abs(PDF(:,2)-quantile));
radius = PDF(nearestIndex,1);
end