color = {[0.7 0 0], [0.043 0.5 0.21], [0.56 0.76 0.24]};
thresh = [0.5 0.5 0.5];
alpha = [0.9 0.99 0.99];
[X,Y,Z]=meshgrid((1:dy)*voxelSizeY,(1:dx)*voxelSizeX,(1:Nz)*voxelSizeZ);
hVol = figure('Name', ['Reconstructed 3D volume: Tpoint_' num2str(tpoint,'%03.f')]);
hVol = plot3DPredictiveRegion(hVol, X, Y, Z, ncRegionMerged, outerCore, innerCore, thresh, color, alpha);

isosurface(X,Y,Z,midPlaneDisp1,0.8);
alpha(0.8)
isosurface(X,Y,Z,midPlaneDisp2,0.8);
alpha(0.8);
axis equal
campos([0,30,25]);