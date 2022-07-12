function [V, R, xg, yg, zg] = EigenVectors3D_FlattenZ(x, y, z, voxelSizeX, voxelSizeY, voxelSizeZ)
%Compute eignen vectors from 3d points represented by x,y and z

xsp=x*voxelSizeX;
ysp=y*voxelSizeY;
zsp=z*voxelSizeZ;
xg=mean(xsp);
yg=mean(ysp);
zg=mean(zsp);
zsp = zg;
Txx=mean(xsp.^2)-xg^2;
Tyy=mean(ysp.^2)-yg^2;
Tzz=mean(zsp.^2)-zg^2;
Txy=mean(xsp.*ysp)-xg*yg;
Txz=mean(xsp.*zsp)-xg*zg;
Tyz=mean(ysp.*zsp)-yg*zg;
Tg=[Txx,Txy,Txz;Txy,Tyy,Tyz;Txz,Tyz,Tzz];

[V,R]=eig(Tg);
end