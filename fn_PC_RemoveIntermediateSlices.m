function [volume] = fn_PC_RemoveIntermediateSlices(intVolume, zFactor)
%Removes interpolates slices from the volume
[dx, dy, dz] = size(intVolume);

newDz = ceil(dz / zFactor);
volume = zeros(dx,dy,newDz);
for i = 1: newDz
    volume(:,:,i) = intVolume(:,:,i*zFactor-zFactor +1);
end
end