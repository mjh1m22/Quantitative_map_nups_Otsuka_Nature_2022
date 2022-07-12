function [inVolume] = threeWayFilling(inVolume, rDisk)
%This functions fills 3D binary volume in three direction xy, yz, xz
[dx, dy, Nz] = size(inVolume);
for zplane = 1: Nz
    inVolume(:,:,zplane)=imdilate(inVolume(:,:,zplane),strel('disk',rDisk,0));
    inVolume(:,:,zplane)=imfill(inVolume(:,:,zplane),'holes');
    for i=1:rDisk
        inVolume(:,:,zplane)=imerode(inVolume(:,:,zplane),strel('diamond',1));
    end
end

for i=1:Nz
    inVolume(:,:,i) = imfill(inVolume(:,:,i),'holes');
end

for i=1:dy
    curSlice = inVolume(:,i,:);
    curSlice = permute(curSlice,[1, 3, 2]);
    curSlice = imfill(curSlice,'holes');
    inVolume(:,i,:) = curSlice;
end

for i=1:dx
    curSlice = inVolume(i,:,:);
    curSlice = permute(curSlice,[2,3,1]);
    curSlice = imfill(curSlice,'holes');
    inVolume(i,:,:) = curSlice;
end
end

