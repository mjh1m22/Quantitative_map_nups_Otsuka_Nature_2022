function [nucRim, centZ] = fn_GetCentralSliceOfRim(nucRim, nucRegion, zFactor)
% Get the cetral plane of nucleus and return nuclear rim this plane and
% the index of the central plan 
[~,~,z] = threeDcoord(find(nucRegion),size(nucRegion,1),size(nucRegion,2));

centZ = round(mean(z));
centZ = round((centZ-1)/zFactor)*zFactor +1;

nucRim(:, :, 1:centZ-1) = 0;
nucRim(:, :, centZ+1:end) = 0;
end