function [ midPlane] = fn_GenerateMidPlane(xg,yg,zg,dx,dy,Nz, voxelSizeX, voxelSizeY, voxelSizeZ, Vma1, Vma2, len)
midPlane = zeros(dx,dy,Nz);
t = -len:0.25:len;
xma1 = t*Vma1(1)+xg;
yma1 = t*Vma1(2)+yg;
zma1 = t*Vma1(3)+zg;

xma2 = t*Vma2(1)+xg;
yma2 = t*Vma2(2)+yg;
zma2 = t*Vma2(3)+zg;

midPlane(:,:,:) = 0;

for i = 1: size(xma1,2)
    yma2cur = round((yma2 - yma1(i) + yg)/voxelSizeY);
    xma2cur = round((xma2 - xma1(i) + xg)/voxelSizeX);
    zma2cur = round((zma2 - zma1(i) + zg)/voxelSizeZ);
    for j = 1: size(xma2cur,2)
        if xma2cur(j)>= 1 && xma2cur(j) <=dx && yma2cur(j)>= 1 && yma2cur(j) <=dy && zma2cur(j)>= 1 && zma2cur(j) <=Nz
            midPlane(xma2cur(j),yma2cur(j),zma2cur(j)) = 1;
        end
    end
end
end

