function [coreProtein] = fn_FindCoreCentroid(nupRegion, protField, protFieldIdx, Vma1, xg,yg,zg, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ, rOrig)
coreProtein = zeros(size(nupRegion));
t = -25:0.05:25;
if sum(sum(sum(protField(:,:,:) == protFieldIdx))) == 0
    return;
end
xma1 = round((t*Vma1(1)+xg)/voxelSizeX);
yma1 = round((t*Vma1(2)+yg)/voxelSizeY);
zma1 = round((t*Vma1(3)+zg)/voxelSizeZ);
coreFlag = 0;
for i = 1: size(xma1,2)
    if xma1(i)>= 1 && xma1(i) <=dx && yma1(i)>= 1 && yma1(i) <=dy && zma1(i)>= 1 && zma1(i) <=Nz
        if (nupRegion(xma1(i), yma1(i),zma1(i)) > 0) && protField(xma1(i), yma1(i),zma1(i)) == protFieldIdx
            xcore = xma1(i);
            ycore = yma1(i);
            zcore = zma1(i);
            coreFlag = 1;
            break;
        end
    end
end
r = round(rOrig);
if coreFlag == 1
    for x=xcore-r:xcore+r
        for y=ycore-r:ycore+r
            for z=zcore-round(r*1.15):zcore+round(r*1.15)
                if x>= 1 && x<=dx && y>= 1 && y<=dy && z>= 1 && z <=Nz && sqrt((x-xcore)^2 +(y-ycore)^2 +(z-zcore)^2) <=rOrig*1.15
                    coreProtein(x,y,z) = 1;
                end
            end
        end 
    end

    coreProtein(protField(:,:,:) ~= protFieldIdx) = 0;
    coreProtein(nupRegion(:,:,:) == 0) = 0;
else
    disp('notFound');
end
end

