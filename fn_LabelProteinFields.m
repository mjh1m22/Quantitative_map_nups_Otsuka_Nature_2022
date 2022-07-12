function [protField1, protField2] = fn_LabelProteinFields(ncRegionCur1, ncRegionCur2, midPlane)
%UNTITLED3 Summary of this function goes here
protField = ones(size(ncRegionCur1));
protField1(:,:,:) = zeros(size(ncRegionCur1));
protField2(:,:,:) = zeros(size(ncRegionCur1));

protField(midPlane(:,:,:) == 1) = 0;

threeDLabel = bwconncomp(protField,26); % Labelling - 3D
numPixels = cellfun(@numel,threeDLabel.PixelIdxList); % Calculate  number of pixels in each connected component
[~,idx] = max(numPixels); % Get the id of biggest connected component
protField1(threeDLabel.PixelIdxList{idx}) = 1;
numPixels(idx) = 0;
[~,idx] = max(numPixels); % Get the id of biggest connected component
protField2(threeDLabel.PixelIdxList{idx}) = 1;

if sum(sum(sum(and(protField2, ncRegionCur2)))) < sum(sum(sum(and(protField1, ncRegionCur2))))
    protField = protField1;
    protField1 = protField2;   
    protField2 = protField;   
end
end

