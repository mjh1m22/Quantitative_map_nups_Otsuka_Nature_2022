function [totSurfArea, protOnSurface] = surfaceArea3D_6Conn_Exclude_Mod(A, E, voxelSizeX, voxelSizeY, voxelSizeZ)
% Construct kernel where we can count the
% number of 6-connected neighbor voxels.
% A is the binary image whose surface is to be estimated
% E is the exlusion mask which is not considered during calculation
conn6Kernel = zeros([3,3,3]);
conn6Kernel(2,2,1) = 1;
conn6Kernel(1,2,2) = 1;
conn6Kernel(2,1,2) = 1;
conn6Kernel(2,3,2) = 1;
conn6Kernel(3,2,2) = 1;
conn6Kernel(2,2,3) = 1;

% For each voxel, determine how many 6-connected neighbors it has.
sumOfFaces = convn(A, conn6Kernel, 'same');

% Find number of exposed faces for each voxel.
surfaceArea = 6 * A - sumOfFaces;
% Mask out zero voxels that have negative exposed faces.
surfaceArea(surfaceArea<0) = 0;
surfaceArea(surfaceArea>1) = 1;
protOnSurface = surfaceArea;
protOnSurface = immultiply(protOnSurface, E);
protOnSurface = imdilate(protOnSurface, ones(3,3,3));
surfaceArea(E>0) = 0;

% Now we simply label the volume and sum up the values of each region.
% binaryVolume = surfaceArea >= 0;
totSurfArea = 4/3 * sum(sum(sum(surfaceArea)))* voxelSizeX * voxelSizeZ;
end