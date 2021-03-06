function stack = tiffread(file_name, indices, varargin)

% tiffread, version 3.1 - October 25, 2012
%
% stack = tiffread;
% stack = tiffread(filename);
% stack = tiffread(filename, indices);
%
% Two options can be specified: 'ReadUnknownTags' and 'DistributeMetaData', eg:
% stack = tiffread(filename, [], 'ReadUnknownTags',1);
% stack = tiffread(filename, [], 'DistributeMetaData', 1);
%
% Reads 8,16,32 bits uncompressed grayscale and (some) color tiff files,
% as well as stacks or multiple tiff images, for example those produced
% by metamorph, Zeiss LSM or NIH-image.
%
% The function can be called with a file name in the current directory,
% or without argument, in which case it pops up a file opening dialog
% to allow for a manual selection of the file.
% If the stacks contains multiples images, reading can be restricted by
% specifying the indices of the desired images (eg. 1:5), or just one index (eg. 2).
%
% The returned value 'stack' is a vector struct containing the images 
% and their meta-data. The length of the vector is the number of images.
% The image pixels values are stored in a field .data, which is a simple
% matrix for gray-scale images, or a cell-array of matrices for color images.
%
% The pixels values are returned in their native (usually integer) format,
% and must be converted to be used in most matlab functions.
%
% Example:
% im = tiffread('spindle.stk');
% imshow( double(im(5).data) );
%
%
%
% Only a fraction of the TIFF standard is supported, but you may extend support
% by modifying this file. If you do so, please return your modification to us,
% such that the added functionality can be redistributed to everyone.
%
% ------------------------------------------------------------------------------
%
% If you would like to acknowledge tiffread2.m in a publication, 
% please cite the article for which the macro was first written:
%
% Dynamic Concentration of Motors in Microtubule Arrays 
% Francois Nedelec, Thomas Surrey and A.C. Maggs
% Physical Review Letters 86: 3192-3195; 2001. 
% DOI:	10.1103/PhysRevLett.86.3192
%
% Thank you!
%
%
% Francois Nedelec
% nedelec -at- embl.de
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
% http://www.embl.org
% http://www.cytosim.org
%
%
% ------------------------------------------------------------------------------
%
% Copyright (C) 1999-2010 Francois Nedelec, 
% with contributions from:
%   Kendra Burbank for the waitbar
%   Hidenao Iwai for the code to read floating point images,
%   Stephen Lang to be more compliant with PlanarConfiguration
%   Jan-Ulrich Kreft for Zeiss LSM support
%   Elias Beauchanp and David Kolin for additional Metamorph support
%   Jean-Pierre Ghobril for requesting that image indices may be specified
%   Urs Utzinger for the better handling of color images, and LSM meta-data
%   O. Scott Sands for support of GeoTIFF tags
%   Benjamin Bratton for Andor tags and more
%   
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details:
% <http://www.gnu.org/licenses/>.
%
% ------------------------------------------------------------------------------
%%

% without argument, we ask the user to choose a file:
if nargin < 1  ||  isempty(file_name)
    [file_name, pathname] = uigetfile('*.tif;*.stk;*.lsm', 'select image file');
    file_name = fullfile(pathname, file_name);
end

if nargin < 2  ||  isempty(indices)
    indices = 1:100000;
end


%% Global Options

parser = inputParser;
parser.addParamValue('ReadUnknownTags',      0);
parser.addParamValue('DistributeMetaData',   1);
parser.addParamValue('ConsolidateStrips',    1);
parser.addParamValue('SkipDissimilarImages', 1);
parser.parse(varargin{:});

opt = parser.Results;


%%

% the structure IMG is returned to the user, while TIF is not.
% so tags usefull to the user should be stored as fields in IMG, while
% those used only internally can be stored in TIF.
% the structure ANDOR has additional header information which is added to
% each plane of the image eventually


global TIF;
TIF = [];
hasAndorHeader = 0;

%% set defaults values :

TIF.SamplesPerPixel  = 1;
TIF.BOS              = 'ieee-le';          %byte order string

if  isempty(strfind(file_name,'.'))
    file_name = [file_name,'.tif'];
end

TIF.file = fopen(file_name,'r','l');
if TIF.file == -1
    stkname = strrep(file_name, '.tif', '.stk');
    TIF.file = fopen(stkname,'r','l');
    if TIF.file == -1
        error(['File "',file_name,'" not found.']);
    else
        file_name = stkname;
    end
end

% obtain the full file path:
[~, file_attrib] = fileattrib(file_name);
file_name = file_attrib.Name;

% obtain the short name:
[~, name, ext] = fileparts(file_name);
image_name = [name, ext];


%% ---- read byte order: II = little endian, MM = big endian

byte_order = fread(TIF.file, 2, '*char');
if ( strcmp(byte_order', 'II') )
    TIF.BOS = 'ieee-le';    % Intel little-endian format
elseif ( strcmp(byte_order','MM') )
    TIF.BOS = 'ieee-be';
else
    error('This is not a TIFF file (no MM or II).');
end


%% ---- read in a number which identifies TIFF format

tiff_id = fread(TIF.file,1,'uint16', TIF.BOS);

if (tiff_id ~= 42)
    error('This is not a TIFF file (missing 42).');
end

%% ---- read the image file directories (IFDs)

ifd_pos   = fread(TIF.file, 1, 'uint32', TIF.BOS);
stack     = [];
hWaitbar  = [];
img_indx  = 0;


while  ifd_pos ~= 0 

    img_indx = img_indx + 1;
    
    clear IMG;
    IMG.file_name = file_name;
    IMG.image_name = image_name;
    IMG.index = img_indx;
 
    % move in the file to the next IFD
    file_seek(ifd_pos);
    %fprintf('reading IFD %i at position %i\n', img_indx, ifd_pos);

    %read in the number of IFD entries
    num_entries = fread(TIF.file,1,'uint16', TIF.BOS);
    %disp(strcat('num_entries =', num2str(num_entries)));

    % store current position:
    entry_pos = ifd_pos+2;
    
    % read the next IFD address:
    file_seek(ifd_pos+12*num_entries+2);
    ifd_pos = fread(TIF.file, 1, 'uint32', TIF.BOS);
 
    % we want to read the first IFD, but not any other one,
    % unless demanded in the indices:
    if img_indx > 1
        if all( img_indx > indices )
            break;
        end
        % skip images that we do not have to read:
        if all( img_indx ~= indices )
            continue;
        end
    end
    
    % read all the IFD entries
    for i = 1:num_entries

        % move to next IFD entry in the file
        file_seek(entry_pos+12*(i-1));

        % read entry tag
        entry_tag = fread(TIF.file, 1, 'uint16', TIF.BOS);
        % read entry
        entry = readIFDentry(entry_tag);
        %disp(strcat('reading entry <',num2str(entry_tag),'>'));


        % not all valid tiff tags have been included, but tags can easily be added to this code
        % See the official list of tags:
        % http://partners.adobe.com/asn/developer/pdfs/tn/TIFF6.pdf
        switch entry_tag
        
            case 254
                TIF.NewSubfiletype = entry.val;
            case 256         % image width = number of column
                IMG.width = entry.val;
            case 257         % image height = number of row
                IMG.height = entry.val;
                TIF.ImageLength = entry.val;
            case 258
                TIF.BitsPerSample = entry.val;
                TIF.BytesPerSample = TIF.BitsPerSample / 8;
                IMG.bits = TIF.BitsPerSample(1);
                %fprintf('BitsPerSample %i %i %i\n', entry.val);
            case 259         % compression
                if ( entry.val ~= 1 )
                    error(['Compression format ', num2str(entry.val),' not supported.']);
                end
            case 262         % photometric interpretation
                TIF.PhotometricInterpretation = entry.val;
                if ( TIF.PhotometricInterpretation == 3 )
                    %warning('tiffread:LookUp', 'Ignoring TIFF look-up
                    %table'); %Commented by Julius
                end
            case 269
                IMG.document_name = entry.val;
            case 270         % general comments:
                IMG.info = entry.val;
            case 271
                IMG.make = entry.val;
            case 273         % strip offset
                TIF.StripOffsets = entry.val;
                TIF.StripNumber = entry.cnt;
                %fprintf('StripNumber = %i, size(StripOffsets) = %i %i\n', TIF.StripNumber, size(TIF.StripOffsets));
            case 274
                % orientation is read, but the matrix is not rotated
                if ( 1 < entry.val ) && ( entry.val < 9 )
                    IMG.orientation = entry.val;
                    keys = {'TopLeft', 'TopRight', 'BottomRight', 'BottomLeft', 'LeftTop', 'RightTop', 'RightBottom', 'LeftBottom'};
                    IMG.orientation_text = keys{entry.val};
                end
            case 277
                TIF.SamplesPerPixel  = entry.val;
                %fprintf('Color image: sample_per_pixel=%i\n',  TIF.SamplesPerPixel);
            case 278
                TIF.RowsPerStrip   = entry.val;
            case 279         % strip byte counts - number of bytes in each strip after any compressio
                TIF.StripByteCounts= entry.val;
            case 282
                IMG.x_resolution   = entry.val;
            case 283
                IMG.y_resolution   = entry.val;
            case 284         %planar configuration describe the order of RGB
                TIF.PlanarConfiguration = entry.val;
            case 296
                IMG.resolution_unit= entry.val;
            case 305
                IMG.software       = entry.val;
            case 306
                IMG.datetime       = entry.val;
            case 315
                IMG.artist         = entry.val;
            case 317        %predictor for compression
                if (entry.val ~= 1); error('unsuported predictor value'); end
            case 320         % color map
                IMG.cmap           = entry.val;
                IMG.colors         = entry.cnt/3;
            case 339
                TIF.SampleFormat   = entry.val;
            
            % ANDOR tags
            case 4864
                hasAndorHeader = 1;
            case 4869       %ANDOR tag: temperature in Celsius when stabilized
               if ~(entry.val == -999)
                   ANDOR.temperature = entry.val;
               end
            case 4876       %exposure time in seconds
                ANDOR.exposureTime   = entry.val;
            case 4878
                ANDOR.kineticCycleTime = entry.val;
            case 4879       %number of accumulations
                ANDOR.nAccumulations = entry.val;
            case 4881
                ANDOR.acquisitionCycleTime = entry.val;
            case 4882       %Readout time in seconds, 1/readoutrate
                ANDOR.readoutTime = entry.val;
            case 4884
                if (entry.val == 9)
                    ANDOR.isPhotonCounting = 1;
                else
                    ANDOR.isPhotonCounting = 0;
                end
            case 4885         %EM DAC level
                ANDOR.emDacLevel = entry.val;
            case 4890
                ANDOR.nFrames = entry.val;
            case 4896
                ANDOR.isFlippedHorizontally = entry.val;
            case 4897
                ANDOR.isFlippedVertically = entry.val;
            case 4898
                ANDOR.isRotatedClockwise = entry.val;
            case 4899
                ANDOR.isRotatedAnticlockwise = entry.val;
            case 4904
                ANDOR.verticalClockVoltageAmplitude = entry.val;
            case 4905
                ANDOR.verticalShiftSpeed = entry.val;
            case 4907
                ANDOR.preAmpSetting = entry.val;
            case 4908         %Camera Serial Number
                ANDOR.serialNumber = entry.val;
            case 4911       %Actual camera temperature when not equal to -999
                if ~(entry.val == -999)
                    ANDOR.unstabilizedTemperature = entry.val;
                end
            case 4912
                ANDOR.isBaselineClamped = entry.val;
            case 4913
                ANDOR.nPrescans = entry.val;
            case 4914
                ANDOR.model = entry.val;
            case 4915
                ANDOR.chipXSize = entry.val;
            case 4916
                ANDOR.chipYSize  = entry.val;
            case 4944
                ANDOR.baselineOffset = entry.val;
            
            case 33550       % GeoTIFF
                IMG.ModelPixelScaleTag = entry.val;
            case 33628       % Metamorph specific data
                IMG.MM_private1 = entry.val;
            case 33629       % this tag identify the image as a Metamorph stack!
                TIF.MM_stack = entry.val;
                TIF.MM_stackCnt = entry.cnt;
            case 33630       % Metamorph stack data: wavelength
                TIF.MM_wavelength = entry.val;
            case 33631       % Metamorph stack data: gain/background?
                TIF.MM_private2 = entry.val;
            
            case 33922       % GeoTIFF
                IMG.ModelTiePointTag = entry.val;
            case 34412       % Zeiss LSM data
                LSM_info = entry.val;
            case 34735       % GeoTIFF
                IMG.GeoKeyDirTag = entry.val;
            case 34737       % GeoTIFF
                IMG.GeoASCII = entry.val;
            case 42113       % GeoTIFF
                IMG.GDAL_NODATA = entry.val;

            case 50838       % Adobe
                IMG.meta_data_byte_counts = entry.val;
            case 50839       % Adobe
                IMG.meta_data = entry.val;
                
            otherwise
                if opt.ReadUnknownTags
                    IMG.(['tag', num2str(entry_tag)])=entry.val;
                    %eval(['IMG.Tag',num2str(entry_tag),'=',entry.val,';']);
                else
                    fprintf( 'Unknown TIFF entry with tag %i (cnt %i)\n', entry_tag, entry.cnt);
                end
        end
        
        % calculate bounding box  if you've got the stuff
        if isfield(IMG, 'ModelPixelScaleTag') && isfield(IMG, 'ModelTiePointTag') && isfield(IMG, 'height')&& isfield(IMG, 'width'),
            IMG.North=IMG.ModelTiePointTag(5)-IMG.ModelPixelScaleTag(2)*IMG.ModelTiePointTag(2);
            IMG.South=IMG.North-IMG.height*IMG.ModelPixelScaleTag(2);
            IMG.West=IMG.ModelTiePointTag(4)+IMG.ModelPixelScaleTag(1)*IMG.ModelTiePointTag(1);
            IMG.East=IMG.West+IMG.width*IMG.ModelPixelScaleTag(1);
        end
        
    end

     
    % do not read image if not demanded
    if all( img_indx ~= indices )
        continue;
    end

    % Planar configuration is not fully supported
    
    %Per tiff spec 6.0 PlanarConfiguration irrelevent if SamplesPerPixel==1
    %Contributed by Stephen Lang
    if (TIF.SamplesPerPixel ~= 1) && ( ~isfield(TIF, 'PlanarConfiguration') || TIF.PlanarConfiguration == 1 )
        error('PlanarConfiguration = 1 is not supported');
    end

    %total number of bytes per image:
    PlaneBytesCnt = IMG.width * IMG.height * TIF.BytesPerSample;

    % try to consolidate the TIFF strips if possible
    
    if opt.ConsolidateStrips
        %Try to consolidate the strips into a single one to speed-up reading:
        BytesCnt = TIF.StripByteCounts(1);

        if BytesCnt < PlaneBytesCnt

            ConsolidateCnt = 1;
            %Count how many Strip are needed to produce a plane
            while TIF.StripOffsets(1) + BytesCnt == TIF.StripOffsets(ConsolidateCnt+1)
                ConsolidateCnt = ConsolidateCnt + 1;
                BytesCnt = BytesCnt + TIF.StripByteCounts(ConsolidateCnt);
                if ( BytesCnt >= PlaneBytesCnt ); break; end
            end

            %Consolidate the Strips
            if ( BytesCnt <= PlaneBytesCnt(1) ) && ( ConsolidateCnt > 1 )
                %fprintf('Consolidating %i stripes out of %i', ConsolidateCnt, TIF.StripNumber);
                TIF.StripByteCounts = [BytesCnt; TIF.StripByteCounts(ConsolidateCnt+1:TIF.StripNumber ) ];
                TIF.StripOffsets = TIF.StripOffsets( [1 , ConsolidateCnt+1:TIF.StripNumber] );
                TIF.StripNumber  = 1 + TIF.StripNumber - ConsolidateCnt;
            end
        end
    end

    
    % Read pixels:
    
    if isfield(TIF, 'MM_stack')

        sel = ( indices <= TIF.MM_stackCnt );
        indices = indices(sel);
        
        if numel(indices) > 1
            hWaitbar = waitbar(0,'Reading images...','Name','TiffRead');
        end

        %this loop reads metamorph stacks:
        for ii = indices

            TIF.StripCnt = 1;
            offset = PlaneBytesCnt * (ii-1);

            %read the image channels
            for c = 1:TIF.SamplesPerPixel
                IMG.data{c} = read_plane(offset, IMG.width, IMG.height, c);
            end

            % print a text timer on the main window, or update the waitbar
            % fprintf('img_read %i img_skip %i\n', img_read, img_skip);
            if ~isempty( hWaitbar )
                waitbar(numel(stack)/numel(indices), hWaitbar);
            end
            
            [ IMG.MM_stack, IMG.MM_wavelength, IMG.MM_private2 ] = splitMetamorph(ii);
            
            IMG.index = ii;
            stack = cat(2, stack, IMG);

        end
        break;

    else
      
        if opt.SkipDissimilarImages && numel(stack) > 1
            if IMG.width ~= stack(1).width || IMG.height ~= stack(1).height
                % setting read_img=0 will skip dissimilar images:
                % comment-out the line below to allow dissimilar stacks
                fprintf('Tiffread skipped %ix%i image\n', width, height);
                read_img = 0;
            end
        end
        
        TIF.StripCnt = 1;
        
        %read the image channels
        for c = 1:TIF.SamplesPerPixel
            IMG.data{c} = read_plane(0, IMG.width, IMG.height, c);
        end
        
        try
            stack = cat(2, stack, IMG);
        catch
            fprintf('Tiffread could not load dissimilar image %i\n', img_indx);
        end

    end
end

%% remove the cell structure if there is always only one channel

flat = 1;
for i = 1:numel(stack)
    if numel(stack(i).data) ~= 1
        flat = 0;
        break;
    end
end

if flat
    for i = 1:numel(stack)
        stack(i).data = stack(i).data{1};
    end
end


%% distribute Andor header to all planes.

if opt.DistributeMetaData  &&  hasAndorHeader

    fieldNames = fieldnames(ANDOR);
    for i = 1:numel(stack)
        for j = 1:size(fieldNames,1)
            stack(i).(fieldNames{j})=ANDOR.(fieldNames{j});
        end
        stack(i).planeNumber = i;
    end
    % set nFrames if it doesn't exist
    if ~ isfield(stack,'nFrames')
        nFrames = numel(stack);
        stack = setfield(stack, {1}, 'nFrames', nFrames);
    end
    
end

%% distribute the MetaMorph info

if opt.DistributeMetaData

    if isfield(TIF, 'MM_stack') && isfield(IMG, 'info') && ~isempty(IMG.info)
        MM = parseMetamorphInfo(IMG.info, TIF.MM_stackCnt);
        for i = 1:numel(stack)
            stack(i).MM = MM(i);
        end
    end
    
    %% duplicate the LSM info
    if exist('LSM_info', 'var')
        for i = 1:numel(stack)
            stack(i).lsm = LSM_info;
        end
    end
    
end

%% return

if ~ exist('stack', 'var')
    stack = [];
end

%clean-up
fclose(TIF.file);
if ~isempty( hWaitbar )
    delete( hWaitbar );
end



end

%%

function file_seek(pos)
global TIF;
status = fseek(TIF.file, pos, -1);
if status == -1
    error('invalid file offset (invalid fseek)');
end
end
        
%% ===========================================================================

function plane = read_plane(offset, width, height, plane_nb)

global TIF;

%return an empty array if the sample format has zero bits
if ( TIF.BitsPerSample(plane_nb) == 0 )
    plane=[];
    return;
end

% fprintf('reading plane %i size %i %i\n', plane_nb, width, height);

% determine the type of data stored in the pixels:
SampleFormat = 1;
if isfield(TIF, 'SampleFormat')  &&  ( plane_nb < length(TIF.SampleFormat) )
    SampleFormat = TIF.SampleFormat(plane_nb);
end

switch( SampleFormat )
    case 1
        classname = sprintf('uint%i', TIF.BitsPerSample(plane_nb));
    case 2
        classname = sprintf('int%i', TIF.BitsPerSample(plane_nb));
    case 3
        if ( TIF.BitsPerSample(plane_nb) == 32 )
            classname = 'single';
        else
            classname = 'double';
        end
    otherwise
        error('unsuported TIFF sample format %i', SampleFormat);
end

% Preallocate a matrix to hold the sample data:
try
    plane = zeros(width, height, classname);
catch
    %compatibility with older matlab versions:
    eval(['plane = ', classname, '(zeros(width, height));']);
end

% Read the strips and concatenate them:
line = 1;
while ( TIF.StripCnt <= TIF.StripNumber )

    strip = read_strip(offset, width, plane_nb, TIF.StripCnt, classname);
    TIF.StripCnt = TIF.StripCnt + 1;

    % copy the strip onto the data
    plane(:, line:(line+size(strip,2)-1)) = strip;

    line = line + size(strip,2);
    if ( line > height )
        break;
    end

end

% Extract valid part of data if needed
if ~all(size(plane) == [width height]),
    plane = plane(1:width, 1:height);
    warning('tiffread:Crop','Cropping data: found more bytes than needed');
end

% transpose the image (otherwise display is rotated in matlab)
plane = plane';

end


%% ================== sub-functions to read a strip ===================

function strip = read_strip(offset, width, plane_nb, stripCnt, classname)

global TIF;

%fprintf('reading strip at position %i\n',TIF.StripOffsets(stripCnt) + offset);
StripLength = TIF.StripByteCounts(stripCnt) ./ TIF.BytesPerSample(plane_nb);

%fprintf( 'reading strip %i\n', stripCnt);
file_seek(TIF.StripOffsets(stripCnt) + offset);

bytes = fread( TIF.file, StripLength, classname, TIF.BOS );

if any( length(bytes) ~= StripLength )
    error('End of file reached unexpectedly.');
end

strip = reshape(bytes, width, StripLength / width);

end


%% ==================sub-functions that reads an IFD entry:===================


function [nbBytes, matlabType] = convertType(tiffType)
switch (tiffType)
    case 1
        nbBytes=1;
        matlabType='uint8';
    case 2
        nbBytes=1;
        matlabType='uchar';
    case 3
        nbBytes=2;
        matlabType='uint16';
    case 4
        nbBytes=4;
        matlabType='uint32';
    case 5
        nbBytes=8;
        matlabType='uint32';
    case 7
        nbBytes=1;
        matlabType='uchar';
    case 11
        nbBytes=4;
        matlabType='float32';
    case 12
        nbBytes=8;
        matlabType='float64';
    otherwise
        error('tiff type %i not supported', tiffType)
end
end

%% ==================sub-functions that reads an IFD entry:===================

function  entry = readIFDentry(entry_tag)

global TIF;
entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.BOS);
entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.BOS);
%disp(['tiffType =', num2str(entry.tiffType),', cnt = ',num2str(entry.cnt)]);

[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32', TIF.BOS);
    %disp(strcat('offset = ', num2str(offset)));
    file_seek(offset);
end


if entry_tag == 33629   % metamorph 'rationals'
    entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.BOS);
elseif entry_tag == 34412  %TIF_CZ_LSMINFO
    entry.val = readLSMinfo;
else
    if entry.tiffType == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);
    else
        entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.BOS);
    end
end

if ( entry.tiffType == 2 );
    entry.val = char(entry.val');
end

end


%% =============distribute the metamorph infos to each frame:
function [MMstack, MMwavelength, MMprivate2] = splitMetamorph(imgCnt)

global TIF;

MMstack = [];
MMwavelength = [];
MMprivate2 = [];

if TIF.MM_stackCnt == 1
    return;
end

left  = imgCnt - 1;

if isfield( TIF, 'MM_stack' )
    S = length(TIF.MM_stack) / TIF.MM_stackCnt;
    MMstack = TIF.MM_stack(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_wavelength' )
    S = length(TIF.MM_wavelength) / TIF.MM_stackCnt;
    MMwavelength = TIF.MM_wavelength(S*left+1:S*left+S);
end

if isfield( TIF, 'MM_private2' )
    S = length(TIF.MM_private2) / TIF.MM_stackCnt;
    MMprivate2 = TIF.MM_private2(S*left+1:S*left+S);
end

end


%% %%  Parse the Metamorph camera info tag into respective fields
% EVBR 2/7/2005, FJN Dec. 2007
function mm = parseMetamorphInfo(info, cnt)

info   = regexprep(info, '\r\n|\o0', '\n');
parse  = textscan(info, '%s %s', 'Delimiter', ':');
tokens = parse{1};
values = parse{2};

first = char(tokens(1,1));

k = 0;
mm = struct('Exposure', zeros(cnt,1));
for i=1:size(tokens,1)
    tok = char(tokens(i,1));
    val = char(values(i,1));
    %fprintf( '"%s" : "%s"\n', tok, val);
    if strcmp(tok, first)
        k = k + 1;
    end
    if strcmp(tok, 'Exposure')
        [v, c, e, pos] = sscanf(val, '%i');
        unit = val(pos:length(val));
        %return the exposure in milli-seconds
        switch( unit )
            case 'ms'
                mm(k).Exposure = v;
            case 's'
                mm(k).Exposure = v * 1000;
            otherwise
                warning('tiffread:Unit', ['Exposure unit "',unit,'" not recognized']);
                mm(k).Exposure = v;
        end
    else
        switch tok
            case 'Binning'
                % Binning: 1 x 1 -> [1 1]
                mm(k).Binning = sscanf(val, '%d x %d')';
            case 'Region'
                mm(k).Region = sscanf(val, '%d x %d, offset at (%d, %d)')';
            otherwise
                field  = regexprep(tok, ' ', '');
                if strcmp(val, 'Off')
                    mm(k).(field) = 0;
                    %eval(['mm(k).',field,'=0;']);
                elseif strcmp(val, 'On')
                    mm(k).(field) = 1;
                    %eval(['mm(k).',field,'=1;']);
                elseif isstrprop(val,'digit')
                    mm(k).(field)= val;
                    %eval(['mm(k).',field,'=str2num(val)'';']);
                else
                    mm(k).(field)= val;
                    %eval(['mm(k).',field,'=val;']);
                end
        end
    end
end

end

%% ==============partial-parse of LSM info:

function R = readLSMinfo()

% Read part of the LSM info table version 2
% this provides only very partial information, since the offset indicate that
% additional data is stored in the file
global TIF;

R.MagicNumber          = sprintf('0x%09X',fread(TIF.file, 1, 'uint32', TIF.BOS));
S.StructureSize        = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionX           = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionY           = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionZ           = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionChannels    = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.DimensionTime        = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.IntensityDataType    = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ThumbnailX           = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.ThumbnailY           = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.VoxelSizeX           = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeY           = fread(TIF.file, 1, 'float64', TIF.BOS);
R.VoxelSizeZ           = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginX              = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginY              = fread(TIF.file, 1, 'float64', TIF.BOS);
R.OriginZ              = fread(TIF.file, 1, 'float64', TIF.BOS);
R.ScanType             = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.SpectralScan         = fread(TIF.file, 1, 'uint16', TIF.BOS);
R.DataType             = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetVectorOverlay  = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetInputLut       = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetOutputLut      = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetChannelColors  = fread(TIF.file, 1, 'uint32', TIF.BOS);
R.TimeInterval         = fread(TIF.file, 1, 'float64', TIF.BOS);
S.OffsetChannelDataTypes = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetScanInformatio = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetKsData         = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetTimeStamps     = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetEventList      = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetRoi            = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetBleachRoi      = fread(TIF.file, 1, 'uint32', TIF.BOS);
S.OffsetNextRecording  = fread(TIF.file, 1, 'uint32', TIF.BOS);

% There are more information stored in this table, which is not read here


%read real acquisition times:
if ( S.OffsetTimeStamps > 0 )
    
    file_seek(S.OffsetTimeStamps);
    
    StructureSize          = fread(TIF.file, 1, 'int32', TIF.BOS);
    NumberTimeStamps       = fread(TIF.file, 1, 'int32', TIF.BOS);
    for i=1:NumberTimeStamps
        R.TimeStamp(i)     = fread(TIF.file, 1, 'float64', TIF.BOS);
    end
    
    %calculate elapsed time from first acquisition:
    R.TimeOffset = R.TimeStamp - R.TimeStamp(1);
    
end


end


