function intVolume= fn_PC_GenIntermediateSlices(volume, zFactor)
%This function generates intermediate slices to convert stack to isotropy
%It uses simple linear interpolation so that neibourhood information is not
%used multiple times.
%It returns interpolated volume (intVolume)
[dx, dy, dz] = size(volume);

newDz = dz * zFactor - zFactor + 1;
intVolume = zeros(dx,dy,newDz);
for i = 1: dz
    intVolume(:,:,i*zFactor-zFactor +1) = volume(:,:,i);
end

for i = 1:dz-1
    for j = 2:zFactor
        intVolume(:,:,i*zFactor-zFactor +j) = intVolume(:,:,i*zFactor-zFactor +1) *(1-(j-1)/zFactor) + intVolume(:,:,i*zFactor +1) * (j-1)/zFactor; 
    end
end