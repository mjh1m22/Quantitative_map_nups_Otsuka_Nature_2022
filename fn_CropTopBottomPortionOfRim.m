function [nucRim] = fn_CropTopBottomPortionOfRim(nucRim, nucRegion, halfPortion)
% Crop top and bottom portion of nucleus and return cropped nuclear rim
activeZ = sum(sum(nucRegion(:,:,:)));
firstZ  = find (activeZ > 0, 1, 'first');
lastZ   = find (activeZ > 0, 1, 'last');
numActiveZ = lastZ - firstZ + 1;
numRemv = round(numActiveZ * halfPortion);

nucRim(:, :, 1 : firstZ+numRemv-1) = 0;
nucRim(:, :, lastZ-numRemv+1 : end) = 0;
end