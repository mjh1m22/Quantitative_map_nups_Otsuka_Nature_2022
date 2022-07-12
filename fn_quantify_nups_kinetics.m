function fn_quantify_nups_kinetics(fp, fn, chr_chan_idx, prot_chan_idx, file_idx,conf)
% This function segments necleoporin channel first. Then it extracts the
% proteins of interest on the surface area of neclear membrane.

% Last update:
% 2016-08-16: Saves segmented nuclear mask both isotropic and non-isotropic
% Detect average intensity from single central (original) plane (xy 3 z 1)
% Number of pixels inside the rim - isotropically (xy pixels, z 3 pixels)

% 2016-07-08 %Includes total intensity and total number of
% pixels in the nuclear rim

% Author: Julius Hossain, EMBL, Heidelberg, Germany: julius.hossain@embl.de
%clear all;
close all; clc
tinit = 1;
if strcmp(fn(end-2:end), 'tif')
    tifFlag = 1;
else
    tifFlag = 0;
end
if tifFlag == 1
    [reader, dim] = getBioformatsReader(fullfile(fp, fn));
    %Get raw stacks containing DNA and POI
    dxOrig = dim.Ny; %Image height
    dyOrig = dim.Nx; % Image width
    NzOrig = dim.Nz; %Number of Z slices in tif files
    %Nc = dim.Nc; %Number of channels
    tfin = dim.Nt;
    %Size of voxels in x, y and z in micro meter
    voxelSizeX = dim.voxelSize(1);
    voxelSizeY = dim.voxelSize(2);
    voxelSizeZ = dim.voxelSize(3);
else
    addpath 'cstruct'
    addpath 'lsm'
    lsminf = lsminfo([fp fn]); %Reads the lsm header
    dxOrig = lsminf.DIMENSIONS(2); %Image height
    dyOrig = lsminf.DIMENSIONS(1); % Image width
    NzOrig = lsminf.DIMENSIONS(3); %Number of Z slices in lsm file
    tfin   = lsminf.DimensionTime;
    voxelSizeX = lsminf.VoxelSizeX * 1e6;
    voxelSizeY = lsminf.VoxelSizeY * 1e6;
    voxelSizeZ = lsminf.VoxelSizeZ * 1e6;
end

bFactor2D = conf.bFactor2D; %Contribution of 2D threshold
minNucVol = conf.minNucVol; % Minimum volume of nucleus
hsize = conf.hsize; %Kernel size for gaussian blur
sigma = conf.sigma; %Sigma for gaussian blur
dSamp = conf.dSamp;
lowRes = conf.lowRes;
numDil = conf.numDil;
numErd = conf.numErd;
in_cor_dm = conf.in_cor_dm;
out_cor_dm = conf.out_cor_dm;
halfPortionChop = conf.halfPortionChop; % how much of nuclear slice will be excluded on the top (30%) and bottom (30%) for rim analyis
refSurfFlag = 0;
midPlaneSE = [0 1 0; 1 1 1; 0 1 0];
protSE = zeros(3,3,3);
protSE(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
protSE(:,:,2) = [0 1 0; 1 1 1; 0 1 0];
protSE(:,:,3) = [0 0 0; 0 1 0; 0 0 0];
%fn = dir([fp '*.lsm']);
idx = find(fp == filesep, 2,'last');
dateStr = date;
outDir = [fp(1:idx(1)) 'Results_' dateStr(8:11) '_' dateStr(4:6) '_' dateStr(1:2) '_' fp(idx(1)+1:idx(2)) fn(1:end-4) filesep];

outDirSegVis = [outDir 'SegVis' filesep];
outDirNucMask = [outDir 'NucMask' filesep];
outDirNucRim = [outDir 'NucRim' filesep];

if ~exist(outDirSegVis)
    mkdir(outDirSegVis);
end
if ~exist(outDirNucMask)
    mkdir(outDirNucMask);
end
if ~exist(outDirNucRim)
    mkdir(outDirNucRim);
end

save(fullfile(outDir, 'configurations.mat'), 'conf');

% Define the list of parameters to be extracted
params_filename = fullfile(outDir, [num2str(file_idx, '%02d') 'extracted_params-prot_channel_' num2str(prot_chan_idx) '.txt']);
paramsName = conf.paramsName;
paramsTable = cell2table(cell(1,length(paramsName)));
paramsTable.Properties.VariableNames = paramsName;

%lsminf = lsminfo([fp fn]); %Reads the lsm header
zinit = 1; %Set higher than one if you want to some of the lower slice
%Size of voxels in x, y and z in micro meter
if lowRes == 1
    voxelSizeX = voxelSizeX*dSamp;
    voxelSizeY = voxelSizeY*dSamp;
    dx = round(dxOrig/dSamp);
    dy = round(dyOrig/dSamp);
else
    dx = dxOrig;
    dy = dyOrig;
end
zFactor = round(voxelSizeZ/voxelSizeX); % zFactor-1 intermediate slice(s) will be generated
Nz = NzOrig * zFactor -zFactor +1;
voxelSizeZ = voxelSizeZ/zFactor; % Update voxelSizeZ based on the number of intermediate slices
voxelSize = voxelSizeX * voxelSizeY * voxelSizeZ; %Volume of a voxel in cubic micro meter
nucStackOrig = zeros(dxOrig,dyOrig,NzOrig-zinit+1); % Allocate memomry for original stack for chromosome
nup107StackOrig = zeros(dxOrig,dyOrig,NzOrig-zinit+1); % Allocate memomry for original stack for chromosome
disp(['Processing: ' fp fn]);

tpoint = tinit;
while tpoint <= tfin
    %Reads input stacks
    tic
    if tifFlag == 1
%         nup107StackOrig = imreadBF([fp fn],zinit:NzOrig,tpoint,1);
%         nucStackOrig    = imreadBF([fp fn],zinit:NzOrig,tpoint,2);
        nucStackOrig  = getTPointBioFormats(reader, dim, tpoint, chr_chan_idx);
        nup107StackOrig = getTPointBioFormats(reader, dim, tpoint, prot_chan_idx);
    else
        for zplane=zinit:NzOrig
            %stacklsm = tiffread31([fp fn], 2*(zplane-zinit+1)-1);
            stacklsm = tiffread31([fp fn],2*((tpoint-1)*NzOrig+zplane-zinit+1)-1);
            nup107StackOrig(:,:,zplane-zinit+1) = stacklsm.data{1};
            nucStackOrig(:,:,zplane-zinit+1) = stacklsm.data{2};
        end
    end
    toc
    
    disp(['Processing: ' num2str(tpoint)]);
    
    if lowRes == 1
        if tpoint == tinit
            nucStackLow = zeros(dx, dy, NzOrig);
            nup107StackLow = zeros(dx, dy, NzOrig);
        end
        for zplane=1:NzOrig
            nucStackLow(:,:,zplane) = imresize(nucStackOrig(:,:,zplane), 1/dSamp, 'bicubic');
            nup107StackLow(:,:,zplane) = imresize(nup107StackOrig(:,:,zplane), 1/dSamp, 'bicubic');
        end
        clear nucStackOrig;
        clear nup107StackOrig
        nucStackOrig = nucStackLow;
        nup107StackOrig = nup107StackLow;
    end
    %Generate isotropic resolution
    nucStack  = fn_PC_GenIntermediateSlices(nucStackOrig, zFactor);
    nup107Stack  = fn_PC_GenIntermediateSlices(nup107StackOrig, zFactor);
    nup107StackBack = nup107Stack;
    clear nucStackOrig;
    clear nup107StackOrig;
    
    if tpoint == tinit
        nucRegion  = zeros(dx,dy,Nz);%Stores the detected chromosome region
        nup107Region  = zeros(dx,dy,Nz);%Stores the detected nup107 region
        nucRegionCur1 = zeros(dx,dy,Nz); %Nucleus region 1
        nucRegionCur2 = zeros(dx,dy,Nz); %Nucleus region 2
        nucRegionPrev1 = zeros(dx,dy,Nz); %Nucleus region in previous frame
        nucRegionSpan1 = zeros(dx,dy,Nz);
        nucRegionSpan2 = zeros(dx,dy,Nz);
        tempVol = zeros(dx,dy,Nz);
        imgCentX = (dx+1)/2;
        imgCentY = (dy+1)/2;
        
        midPlane1 = zeros(dx,dy,Nz);
        midPlane2 = zeros(dx,dy,Nz);
        midPlaneDisp1 = zeros(dx,dy,Nz);
        midPlaneDisp2 = zeros(dx,dy,Nz);
        midPlane3 = zeros(dx,dy,Nz);
        mergedProtOnSurface = zeros(dx,dy,Nz);
        outerCore = zeros(dx,dy,Nz);
        innerCore = zeros(dx,dy,Nz);
        forcedNuc1 = 0;
        forcedNuc2 = 0;
        nucRimDisp = zeros(dx,dy);
    end
    
    nucStack = imgaussian(nucStack, sigma, hsize);
    nucStackNorm = nucStack;
    normFlag = 1;
    
    while (normFlag == 1)
        [nucThresh3D, nucThresh2D] = fn_PC_GetThresholds(nucStackNorm);
        for i = 1:Nz
            nucRegion(:,:,i) = double(nucStack(:,:,i) >= (nucThresh3D * (1-bFactor2D) + nucThresh2D(i) * bFactor2D));
        end
        
        threeDLabel = bwconncomp(nucRegion,18); % Labelling - 3D
        numPixels = cellfun(@numel,threeDLabel.PixelIdxList); % Calculate  number of pixels in each connected component
        [biggestPix,idx] = max(numPixels);
        biggestVol = biggestPix * voxelSize;
        if biggestVol < minNucVol
            nucStackNorm(threeDLabel.PixelIdxList{idx}) = nucThresh3D;
        else
            normFlag = 0;
        end
    end
    
    threeDLabel = bwconncomp(nucRegion,18); % Labelling - 3D
    numPixels = cellfun(@numel,threeDLabel.PixelIdxList); % Calculate  number of pixels in each connected component
    
    for i = 1 : numel(numPixels)
        if numPixels(i) * voxelSize < minNucVol*0.66
            nucRegion(threeDLabel.PixelIdxList{i}) = 0;
        end
    end
    
    threeDLabel = bwconncomp(nucRegion,18); % Labelling - 3D
    numPixels = cellfun(@numel,threeDLabel.PixelIdxList); % Calculate  number of pixels in each connected component
    
    %%
    if tpoint == tinit
        %Use location information to select the nucleus of interest
        stat = regionprops(threeDLabel,'Centroid');
        %Normalize the number of pixels based on three properties
        for i = 1: numel(stat)
            bDist = ((imgCentX-stat(i).Centroid(2))^2 + (imgCentY - stat(i).Centroid(1))^2)^2; %Special case you expect nc at centre;
            numPixels(i) = numPixels(i) * 1/bDist;
        end
        %Select the best - discard the one with artificial chromosome
        [~,idx] = max(numPixels);
        nucRegionCur1(threeDLabel.PixelIdxList{idx}) = 1;
        curCentX1 = stat(idx).Centroid(2);
        curCentY1 = stat(idx).Centroid(1);
        numNuc = 1;
    else
        %Use location information to select the nucleus of interest
        stat = regionprops(threeDLabel,'Centroid');
        %Normalize the number of pixels based on three properties
        for i = 1: numel(stat)
            bDist = (((prevCentX1-stat(i).Centroid(2))^2 + (prevCentY1 - stat(i).Centroid(1))^2))^3;
            numPixels(i) = numPixels(i) * 1/bDist;
        end
        %Select the best - discard the one with artificial chromosome
        nucRegionCur1(:,:,:) = 0;
        [~,idx] = max(numPixels);
        nucRegionCur1(threeDLabel.PixelIdxList{idx}) = 1;
        numPixels(idx) = 0;
        %         volNcCur1 = sum(sum(sum(nucRegionCur1)));
        %         if volNcCur1 > sum(sum(sum(nucRegionPrev1)))*2.5 || volNcCur1 < sum(sum(sum(nucRegionPrev1)))/3
        %             continue;
        %         end
        
        curCentX1 = stat(idx).Centroid(2);
        curCentY1 = stat(idx).Centroid(1);
        nuc1Idx = idx;
        
        if numNuc == 1
            if sum(sum(sum(nucRegionCur1))) < sum(sum(sum(nucRegionPrev1)))*0.66
                tempVol(:,:,:) = 0;
                [~,idx] = max(numPixels);
                if idx ~= nuc1Idx
                    tempVol(threeDLabel.PixelIdxList{idx}) = 1;
                    
                    distPrevCh1 = sqrt((prevCentX1-curCentX1)^2 + (prevCentY1-curCentY1)^2);
                    distPrevCh2 = sqrt((prevCentX1-stat(idx).Centroid(2))^2 + (prevCentY1-stat(idx).Centroid(1))^2);
                    vol1 = sum(sum(sum(nucRegionCur1)));
                    vol2 = sum(sum(sum(tempVol)));
                    totalVol = vol1 + vol2;
                    if totalVol <= sum(sum(sum(nucRegionPrev1)))*1.3 && distPrevCh1 <= distPrevCh2*4 && distPrevCh2 <= distPrevCh1*4 && vol1 >= vol2*0.6 && vol2 >= vol1*0.6
                        nucRegionCur2(threeDLabel.PixelIdxList{idx}) = 1;
                        curCentX2 = stat(idx).Centroid(2);
                        curCentY2 = stat(idx).Centroid(1);
                        curCentZ2 = stat(idx).Centroid(3);
                        numNuc = 2;
                    end
                end
            end
        else
            if sum(sum(sum(nucRegionCur1))) <= sum(sum(sum(nucRegionPrev1)))*1.9
                numPixels = cellfun(@numel,threeDLabel.PixelIdxList);
                for i = 1: numel(stat)
                    if i == nuc1Idx
                        numPixels(i) = 0;
                    end
                    bDist = (((prevCentX2-stat(i).Centroid(2))^2 + (prevCentY2 - stat(i).Centroid(1))^2))^2;
                    numPixels(i) = numPixels(i) * 1/bDist;
                end
                %Select the best - discard the one with artificial chromosome
                nucRegionCur2(:,:,:) = 0;
                [~,idx] = max(numPixels);
                nucRegionCur2(threeDLabel.PixelIdxList{idx}) = 1;
                %numPixels(idx) = 0;
                
                curCentX2 = stat(idx).Centroid(2);
                curCentY2 = stat(idx).Centroid(1);
                curCentZ2 = stat(idx).Centroid(3);
            end
        end
    end
    
    prevCentX1 = curCentX1;
    prevCentY1 = curCentY1;
    nucRegionPrev1 = nucRegionCur1;
    
    if numNuc == 2
        prevCentX2 = curCentX2;
        prevCentY2 = curCentY2;
    end
    
    nucRegionCur1 = imclose(nucRegionCur1, protSE);
    for i=1:Nz
        nucRegionCur1(:,:,i) = imfill(nucRegionCur1(:,:,i), 'holes');
    end
    nucRegionCur1 = fn_SmoothBW3D(nucRegionCur1, 7);
    
    %     if tpoint > tinit
    %         nucRegionCur1(nucRegionSpan(:,:,:) == 0) = 0;
    %     end
    
    totPixel1 = sum(sum(sum(nucRegionCur1)));
    volume1 = totPixel1 * voxelSize;
    if volume1== 0
        nucRegionCur1 = nucRegionPrev1;
    end
    nucRegionPrev1 = nucRegionCur1;
    totSurfArea1 = surfaceArea3D_6Conn(nucRegionCur1, voxelSizeX, voxelSizeY, voxelSizeZ);
    try
        totSurfArea1Soft = imSurface_Modified(nucRegionCur1, voxelSizeX, voxelSizeY, voxelSizeZ);
    catch
        totSurfArea1Soft = totSurfArea1;
        disp('Could not calculate Surface Area');
    end
    %end
    if numNuc == 2
        nucRegionCur2 = imclose(nucRegionCur2, protSE);
        for i=1:Nz
            nucRegionCur2(:,:,i) = imfill(nucRegionCur2(:,:,i), 'holes');
        end
        nucRegionCur2 = fn_SmoothBW3D(nucRegionCur2, 7);
%         if tpoint > tinit
%             nucRegionCur2(nucRegionSpan(:,:,:) == 0) = 0;
%         end
        
        totPixel2 = sum(sum(sum(nucRegionCur2)));
        volume2 = totPixel2 * voxelSize;
        if volume2== 0
            nucRegionCur2 = nucRegionPrev2;
        end
        totSurfArea2 = surfaceArea3D_6Conn(nucRegionCur2, voxelSizeX, voxelSizeY, voxelSizeZ);
        try
            totSurfArea2Soft = imSurface_Modified(nucRegionCur2, voxelSizeX, voxelSizeY, voxelSizeZ);
        catch
            totSurfArea2Soft = totSurfArea2;
            disp('Could not calculate Surface Area');
        end
        %% calculate total protein on the rim of nuc2
        %protRegionNuc2 = nucRegionCur2;
        %protRegionNuc2 = double(imdilate(protRegionNuc2,protSE) - imerode(imerode(protRegionNuc2, protSE), protSE));
        protRegionNuc2 = fn_generate_nuclear_rim(nucRegionCur2, numDil, numErd, protSE);
        protStack = nup107Stack;
        protStack(protRegionNuc2(:,:,:) == 0) = 0;
        %protStackNonIso  = fn_PC_RemoveIntermediateSlices(protStack, zFactor);
        %protRegionNuc2NonIso = fn_PC_RemoveIntermediateSlices(protRegionNuc2, zFactor);
        %         totIntNuc2 = sum(sum(sum(protStack(protRegionNuc2(:,:,:)>0))));
        %         totNumPixNuc2 = sum(sum(sum(protRegionNuc2(:,:,:)>0)));
        protRegionNuc2ms = fn_CropTopBottomPortionOfRim(protRegionNuc2, nucRegionCur2, halfPortionChop);
        [avgIntNuc2ms, totPixNuc2ms] = fn_calculate_rim_intensity(protRegionNuc2ms, protStack);
    end
    %% calculate total protein on the rim of nuc1
    %protRegionNuc1 = nucRegionCur1;
    %protRegionNuc1 = double(imdilate(protRegionNuc1,protSE) - imerode(imerode(protRegionNuc1, protSE), protSE));
    protRegionNuc1 = fn_generate_nuclear_rim(nucRegionCur1, numDil, numErd, protSE);
    protStack = nup107Stack;
    protStack(protRegionNuc1(:,:,:) == 0) = 0;
    
    protRegionNuc1ms = fn_CropTopBottomPortionOfRim(protRegionNuc1, nucRegionCur1, halfPortionChop);
    [avgIntNuc1ms, totPixNuc1ms] = fn_calculate_rim_intensity(protRegionNuc1ms, protStack);
    %     totIntNuc1 = sum(sum(sum(protStack(protRegionNuc1(:,:,:)>0))));
    %     totNumPixNuc1 = sum(sum(sum(protRegionNuc1(:,:,:)>0)));
    
    %%
    nucRegionPrev2 = nucRegionCur2;
    prevDet = numNuc;
    
    nucRegionSpan = imdilate(double(or(nucRegionCur1, nucRegionCur2)), ones(9,9,5));
    
    % Process nup3582alpha
    extNucRegion1 = imdilate(nucRegionCur1, ones(5,5,3));
    extNucRegion2 = imdilate(nucRegionCur2, ones(5,5,3));
    %nup358Stack(extNucRegion1(:,:,:) == 0 & extNucRegion2(:,:,:) == 0) = 0;
    
    %protRegion = double(nup358Stack(:,:,:) >= refThresh);
    nucRegion = double(nucRegionCur1 | nucRegionCur2);
    %nup107Region1 = nucRegionCur1;
    nup107Region1 = fn_generate_nuclear_rim(nucRegionCur1, numDil, numErd, protSE);
    nup107Region= nup107Region1;
    %nup107Region1 = double(imdilate(nup107Region1,protSE) - imerode(imerode(nup107Region1, protSE), protSE));
    mergedProtOnSurface(:,:,:) = 0;
    nucRimDisp(:,:) = 0;
    if numNuc ==2
        [x,y,z] = threeDcoord(find(nucRegionCur1),dx,dy);
        [V1, R1, xg1, yg1, zg1] = EigenVectors3D_FlattenZ(x, y, z, voxelSizeX, voxelSizeY, voxelSizeZ);
        %         if (abs(V1(3,3)) < 0.7 && abs(V1(3,2)) < 0.7) || forcedNuc1 == 1
        %             xlswrite(xlsFileName,V1(3,2),1,strcat('N', num2str(tpoint+1)));
        %             V1(:,2) = [0 0 1];
        %             %forcedNuc1 = 1;
        %         end
        %         V1(:,1) = cross(V1(:,3), V1(:,2));
        midPlane1     = fn_GenerateMidPlane(xg1,yg1,zg1,dx,dy,Nz, voxelSizeX, voxelSizeY, voxelSizeZ, V1(:,3), V1(:,1), round(50/voxelSizeX));
        midPlaneDisp1 = fn_GenerateMidPlane(xg1,yg1,zg1,dx,dy,Nz, voxelSizeX, voxelSizeY, voxelSizeZ, V1(:,3), V1(:,1), round(2/voxelSizeX));
        midPlane1 = imdilate(midPlane1, ones(3,3,3));
        midPlaneDisp1 = imdilate(midPlaneDisp1, midPlaneSE);
        
        [x,y,z] = threeDcoord(find(nucRegionCur2),dx,dy);
        [V2, R2, xg2, yg2, zg2] = EigenVectors3D_FlattenZ(x, y, z, voxelSizeX, voxelSizeY, voxelSizeZ);
        %Vma2 = V2(:,1);
        %         if (abs(V2(3,3)) < 0.7 && abs(V2(3,2)) < 0.7) || forcedNuc2 == 1
        %             xlswrite(xlsFileName,V2(3,2),1,strcat('O', num2str(tpoint+1)));
        %             V2(:,2) = [0 0 1];
        %             %forcedNuc2 = 1;
        %         end
        %         V2(:,1) = cross(V2(:,3), V2(:,2));
        
        midPlane2     = fn_GenerateMidPlane(xg2,yg2,zg2,dx,dy,Nz, voxelSizeX, voxelSizeY, voxelSizeZ, V2(:,3), V2(:,1), round(50/voxelSizeX));
        midPlaneDisp2 = fn_GenerateMidPlane(xg2,yg2,zg2,dx,dy,Nz, voxelSizeX, voxelSizeY, voxelSizeZ, V2(:,3), V2(:,1), round(2/voxelSizeX));
        midPlane2 = imdilate(midPlane2, ones(3,3,3));
        midPlaneDisp2 = imdilate(midPlaneDisp2, midPlaneSE);
        
        if refSurfFlag == 0
            [protField1, protField2] = fn_LabelProteinFields(nucRegionCur1, nucRegionCur2, midPlane1);
            [protField4, protField3] = fn_LabelProteinFields(nucRegionCur2, nucRegionCur1, midPlane2);
        else
            [protField1, protField2] = fn_LabelProteinFields(protField1, protField2, midPlane1);
            [protField3, protField4] = fn_LabelProteinFields(protField3, protField4, midPlane2);
        end
        
        if refSurfFlag == 0
            refSurfArea1 = totSurfArea1Soft;
            refSurfArea2 = totSurfArea2Soft;
            refSurfFlag = 1;
        end
        
        %nup107Region = nucRegion;
        %nup107Region2 = nucRegionCur2;
        nup107Region = fn_generate_nuclear_rim(nucRegion, numDil, numErd, protSE);
        nup107Region2 = fn_generate_nuclear_rim(nucRegionCur2, numDil, numErd, protSE);
        %nup107Region = double(imdilate(nup107Region,protSE) - imerode(imerode(nup107Region, protSE), protSE));
        %nup107Region2 = double(imdilate(nup107Region2,protSE) - imerode(imerode(nup107Region2, protSE), protSE));
        
        rFactor1 = sqrt(totSurfArea1Soft/refSurfArea1);
        rFactor2 = sqrt(totSurfArea2Soft/refSurfArea2);
        outerCore1 = fn_FindCoreCentroid(nup107Region1, protField1, 1, V1(:,2), xg1,yg1,zg1, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(out_cor_dm/voxelSizeX)*rFactor1);
        outerCore4 = fn_FindCoreCentroid(nup107Region2, protField4, 1, V2(:,2), xg2,yg2,zg2, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(out_cor_dm/voxelSizeX)*rFactor2);
        outerCore = double(or(outerCore1, outerCore4));
        
        innerCore2 = fn_FindCoreCentroid(nup107Region1, protField2, 1, V1(:,2), xg1,yg1,zg1, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(in_cor_dm/voxelSizeX)*rFactor1);
        innerCore3 = fn_FindCoreCentroid(nup107Region2, protField3, 1, V2(:,2), xg2,yg2,zg2, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(in_cor_dm/voxelSizeX)*rFactor2);
        innerCore = double(or(innerCore2, innerCore3));
        
        nup107Stack(nup107Region(:,:,:) == 0) = 0;
        
        outerCore1ms = fn_CropTopBottomPortionOfRim(outerCore1, nucRegionCur1, halfPortionChop);
        [outerCore1ss, ~] = fn_GetCentralSliceOfRim(outerCore1, nucRegionCur1, zFactor);
        %         totOuterInt = sum(sum(sum(nup107Stack(outerCore1(:,:,:)>0))));
        %         totPixOuterCore1ms = sum(sum(sum(outerCore1(:,:,:)>0)));
        %         avgIntOuterCore1ms = totOuterInt/totPixOuterCore1ms;
        [avgIntOuterCore1ms, totPixOuterCore1ms] = fn_calculate_rim_intensity(outerCore1ms, nup107Stack);
        [avgIntOuterCore1ss, totPixOuterCore1ss] = fn_calculate_rim_intensity(outerCore1ss, nup107Stack);
        
        
        %         totInnerInt = sum(sum(sum(nup107Stack(innerCore2(:,:,:)>0))));
        %         totInnerPix107_2 = sum(sum(sum(innerCore2(:,:,:)>0)));
        %         innerAvgInt107_2 = totInnerInt/totInnerPix107_2;
        innerCore2ms = fn_CropTopBottomPortionOfRim(innerCore2, nucRegionCur1, halfPortionChop);
        [innerCore2ss, ~] = fn_GetCentralSliceOfRim(innerCore2, nucRegionCur1, zFactor);
        [avgIntInnerCore2ms, totPixInnerCore2ms] = fn_calculate_rim_intensity(innerCore2ms, nup107Stack);
        [avgIntInnerCore2ss, totPixInnerCore2ss] = fn_calculate_rim_intensity(innerCore2ss, nup107Stack);
        
        %         totOuterInt = sum(sum(sum(nup107Stack(outerCore4(:,:,:)>0))));
        %         totOuterPix107_4 = sum(sum(sum(outerCore4(:,:,:)>0)));
        %         outerAvgInt107_4 = totOuterInt/totOuterPix107_4;
        outerCore4ms = fn_CropTopBottomPortionOfRim(outerCore4, nucRegionCur2, halfPortionChop);
        [outerCore4ss, ~] = fn_GetCentralSliceOfRim(outerCore4, nucRegionCur2, zFactor);
        [avgIntOuterCore4ms, totPixOuterCore4ms] = fn_calculate_rim_intensity(outerCore4ms, nup107Stack);
        [avgIntOuterCore4ss, totPixOuterCore4ss] = fn_calculate_rim_intensity(outerCore4ss, nup107Stack);
        
        %         totInnerInt = sum(sum(sum(nup107Stack(innerCore3(:,:,:)>0))));
        %         totInnerPix107_3 = sum(sum(sum(innerCore3(:,:,:)>0)));
        %         innerAvgInt107_3 = totInnerInt/totInnerPix107_3;
        innerCore3ms = fn_CropTopBottomPortionOfRim(innerCore3, nucRegionCur2, halfPortionChop);
        [innerCore3ss, ~] = fn_GetCentralSliceOfRim(innerCore3, nucRegionCur2, zFactor);
        [avgIntInnerCore3ms, totPixInnerCore3ms] = fn_calculate_rim_intensity(innerCore3ms, nup107Stack);
        [avgIntInnerCore3ss, totPixInnerCore3ss] = fn_calculate_rim_intensity(innerCore3ss, nup107Stack);
        
        %For display
        innerCore = imdilate(innerCore, protSE);
        outerCore = imdilate(outerCore, protSE);
        
        %% Calculate non core intensity
        outerCoreEx1 = fn_FindCoreCentroid(nup107Region1, protField1, 1, V1(:,2), xg1,yg1,zg1, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(out_cor_dm*2/voxelSizeX)*rFactor1);
        outerCoreEx2 = fn_FindCoreCentroid(nup107Region2, protField4, 1, V2(:,2), xg2,yg2,zg2, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(out_cor_dm*2/voxelSizeX)*rFactor2);
        outerCoreEx = double(or(outerCoreEx1, outerCoreEx2));
        
        innerCoreEx1 = fn_FindCoreCentroid(nup107Region1, protField2, 1, V1(:,2), xg1,yg1,zg1, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(in_cor_dm*2/voxelSizeX-1)*rFactor1);
        innerCoreEx2 = fn_FindCoreCentroid(nup107Region2, protField3, 1, V2(:,2), xg2,yg2,zg2, dx,dy, Nz, voxelSizeX, voxelSizeY, voxelSizeZ,round(in_cor_dm*2/voxelSizeX-1)*rFactor2);
        innerCoreEx = double(or(innerCoreEx1, innerCoreEx2));
        nonCore12 = double(or (outerCoreEx1, innerCoreEx1));
        nonCore12 = nup107Region1 - nonCore12; %Corrected 2016-07-11: it was nup107Region instead
        nonCore12(nonCore12<0) = 0;
        
        nonCore34 = double(or (outerCoreEx2, innerCoreEx2));
        nonCore34 = nup107Region2 - nonCore34; %Corrected 2016-07-11: it was nup107Region instead
        nonCore34(nonCore34<0) = 0;
        %nonCore = imdilate(nonCore, protSE);
        nonCore12ms = fn_CropTopBottomPortionOfRim(nonCore12, nucRegionCur1, halfPortionChop);
        [nonCore12ss, ~] = fn_GetCentralSliceOfRim(nonCore12, nucRegionCur1, zFactor);
        [avgIntNonCore12ms, totPixNonCore12ms] = fn_calculate_rim_intensity(nonCore12ms, nup107Stack);
        [avgIntNonCore12ss, totPixNonCore12ss] = fn_calculate_rim_intensity(nonCore12ss, nup107Stack);
        %         totNonInt = sum(sum(sum(nup107Stack(nonCore12(:,:,:)>0))));
        %         totNonPix107_12 = sum(sum(sum(nonCore12(:,:,:)>0)));
        %         nonAvgInt107_12 = totNonInt/totNonPix107_12;
        
        nonCore34ms = fn_CropTopBottomPortionOfRim(nonCore34, nucRegionCur2, halfPortionChop);
        [nonCore34ss, ~] = fn_GetCentralSliceOfRim(nonCore34, nucRegionCur2, zFactor);
        [avgIntNonCore34ms, totPixNonCore34ms] = fn_calculate_rim_intensity(nonCore34ms, nup107Stack);
        [avgIntNonCore34ss, totPixNonCore34ss] = fn_calculate_rim_intensity(nonCore34ss, nup107Stack);
        %         totNonInt = sum(sum(sum(nup107Stack(nonCore34(:,:,:)>0))));
        %         totNonPix107_34 = sum(sum(sum(nonCore34(:,:,:)>0)));
        %         nonAvgInt107_34 = totNonInt/totNonPix107_34;
        
        %% Calculating average intensity of nup107 in the central slice of nuc2
        [nucRim2ss, centZ2] = fn_GetCentralSliceOfRim(nup107Region2, nucRegionCur2, zFactor);
        [avgIntNuc2ss, ~] = fn_calculate_rim_intensity(nucRim2ss, nup107Stack);
        %         totIntNuc2 = sum(sum(sum(nup107Stack(nucRim2ss(:,:,:)>0))));
        %         totPixNuc2 = sum(sum(sum(nucRim2ss(:,:,:)>0)));
        %         avgIntNuc2ss = totIntNuc2/totPixNuc2;
        isonucRim2ss = nucRegionCur2;
        isonucRim2ss(:,:,1:centZ2-2) = 0;
        isonucRim2ss(:,:,centZ2+2:end) = 0;
        totPixNuc2ss = sum(sum(sum(isonucRim2ss(:,:,:)>0)));
        
        nucRimDisp2 = nup107Stack(:,:,centZ2);
        nucRimDisp2(nucRim2ss(:,:,centZ2) == 0) = 0;
    end
    pix_ratio_lhi_vs_uhi =0;
    %[pix_ratio_lhi_vs_uhi] = fn_pix_ratio_lower_half_int_vs_upper_half_int(nup107Region, nup107Stack, 0.5);
    %% Calculating average intensity of nup107 in the central slice of nuc1
    [nucRim1ss, centZ1] = fn_GetCentralSliceOfRim(nup107Region1, nucRegionCur1, zFactor);
    [avgIntNuc1ss, ~] = fn_calculate_rim_intensity(nucRim1ss, nup107Stack);
    %     totIntNuc1ss = sum(sum(sum(nup107Stack(nucRim1(:,:,:)>0))));
    %     totPixNuc1ss = sum(sum(sum(nucRim1(:,:,:)>0)));
    %     avgIntNuc1ss = totIntNuc1ss/totPixNuc1ss;
    isoNucRim1ss = nucRegionCur1;
    isoNucRim1ss(:,:,1:centZ1-2) = 0;
    isoNucRim1ss(:,:,centZ1+2:end) = 0;
    totPixNuc1ss = sum(sum(sum(isoNucRim1ss(:,:,:)>0)));
    
    nucRimDisp1 = nup107Stack(:,:,centZ1);
    nucRimDisp1(nucRim1ss(:,:,centZ1) == 0) = 0;
    if numNuc == 1
        nucRimDisp = nucRimDisp1;
    else
        nucRimDisp = imadd(nucRimDisp1, nucRimDisp2);
    end
    
    %% Display 3D reconstructed volumes and save files
    nucRegionCur2(1,1,1) = 1; nucRegionCur2(1,end,1) = 1; nucRegionCur2(1,1,end) = 1;
    nucRegionCur2(end,1,1) = 1; nucRegionCur2(end,1,end) = 1; nucRegionCur2(1,end,end) = 1;
    nucRegionCur2(end,end,end) = 1; nucRegionCur2(end,end,1) = 1;
    ncRegionMerged = double(nucRegionCur1 | nucRegionCur2);
    %     midPlaneSpan = imdilate(ncRegionMerged, ones(9,9,9));
    %     midPlaneDisp1(midPlaneSpan(:,:,:) == 0) = 0;
    %     midPlaneDisp2(midPlaneSpan(:,:,:) == 0) = 0;
    
    [X,Y,Z]=meshgrid((1:dy)*voxelSizeY,(1:dx)*voxelSizeX,(1:Nz)*voxelSizeZ);
    hVol = figure('Name', ['Reconstructed 3D volume: Tpoint_' num2str(tpoint,'%03.f')]);
    isosurface(X,Y,Z,ncRegionMerged,0.9);
    alpha(0.5)
    isosurface(X,Y,Z,outerCore,0.7);
    alpha(.7)
    isosurface(X,Y,Z,innerCore,0.6);
    alpha(.6)
    isosurface(X,Y,Z,midPlaneDisp1,0.8);
    alpha(0.8)
    isosurface(X,Y,Z,midPlaneDisp2,0.8);
    alpha(0.8);
    axis equal
    campos([0,30,25]);
    
    saveFilename = [outDirSegVis fn(1:end-4) num2str(tpoint,'%03.f') '.tif'];
    saveas(hVol, saveFilename, 'tif');
    delete(hVol);
    
    nucMaskIso = nucRegionCur1;
    nucMaskIso(nucRegionCur2(:,:,:)>0) = 2;
    nucMaskNonIso = fn_PC_RemoveIntermediateSlices(nucMaskIso, zFactor);
    saveFilename = [outDirNucMask fn(1:end-4) num2str(tpoint,'%03.f') '.mat'];
    save(saveFilename, 'numNuc', 'nucMaskIso', 'nucMaskNonIso');
    
    saveFilename = [outDirNucRim fn(1:end-4) num2str(tpoint,'%03.f') '.tif'];
    imwrite(uint16(nucRimDisp), saveFilename, 'tif');
    
    saveFilename = [outDirNucRim fn(1:end-4) num2str(tpoint,'%03.f') '_1' '.tif'];
    imwrite(uint16(nup107StackBack(:,:,centZ1)), saveFilename, 'tif');
    if numNuc == 2
        saveFilename = [outDirNucRim fn(1:end-4) num2str(tpoint,'%03.f') '_2' '.tif'];
        imwrite(uint16(nup107StackBack(:,:,centZ2)), saveFilename, 'tif');
    end

    if numNuc == 2
        newTable = cell2table({tpoint, volume1, totSurfArea1Soft, volume2, totSurfArea2Soft, avgIntOuterCore1ms, avgIntInnerCore2ms, avgIntInnerCore3ms, avgIntOuterCore4ms, avgIntNonCore12ms, avgIntNonCore34ms,...
                   totPixOuterCore1ms, totPixInnerCore2ms, totPixInnerCore3ms, totPixOuterCore4ms, totPixNonCore12ms, totPixNonCore34ms,  avgIntNuc1ms, totPixNuc1ms, avgIntNuc2ms, totPixNuc2ms,...
                   avgIntOuterCore1ss, avgIntInnerCore2ss, avgIntInnerCore3ss, avgIntOuterCore4ss, avgIntNonCore12ss, avgIntNonCore34ss, totPixOuterCore1ss, totPixInnerCore2ss, totPixInnerCore3ss,...
                   totPixOuterCore4ss, totPixNonCore12ss, totPixNonCore34ss, avgIntNuc1ss, totPixNuc1ss, avgIntNuc2ss, totPixNuc2ss, pix_ratio_lhi_vs_uhi});
    else
        newTable = cell2table({tpoint, volume1, totSurfArea1Soft, 0, 0, 0, 0, 0, 0, 0,0,...
                   0, 0, 0, 0, 0, 0,  avgIntNuc1ms, totPixNuc1ms, 0, 0,...
                   0, 0, 0, 0, 0, 0, 0, 0, 0,...
                   0, 0, 0, avgIntNuc1ss, totPixNuc1ss, 0, 0, pix_ratio_lhi_vs_uhi});
    end
    newTable.Properties.VariableNames = paramsName;
    paramsTable = cat(1, paramsTable, newTable);
    
    if tpoint == tinit
        paramsTable(1,:) = [];
    end
    writetable(paramsTable, params_filename, 'Delimiter','\t');
    close all;
    tpoint = tpoint + 1;
end
clear all;
end