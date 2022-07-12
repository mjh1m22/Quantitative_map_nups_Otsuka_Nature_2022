function [avgInt, totPix] = fn_calculate_rim_intensity(nucRim, nupStack)
%UNTITLED3 Summary of this function goes here
totInt = sum(sum(sum(nupStack(nucRim(:,:,:)>0))));
totPix = sum(sum(sum(nucRim(:,:,:)>0)));
avgInt = totInt/totPix;
end

