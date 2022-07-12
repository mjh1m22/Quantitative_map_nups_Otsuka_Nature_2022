function [conf] = fn_configurations()
%% This fuction sets the running options and parameters
% conf: structure returning running options and parameters

% Set the input directory containing tif files to process
conf.in_dir = 'C:\Research_materials\Data\NPC_core_non_core';
%Channel index of chromatin
conf.chr_chan_idx = 2;
%Channel index of nup
conf.prot_chan_idx = 1;

conf.bFactor2D = 0.15; % Contribution of 2D threshold
conf.minNucVol = 175; % Minimum volume of nucleus in cubic micrometer
conf.hsize = 3; % Kernel size for gaussian blur
conf.sigma = 2; % Sigma for gaussian blur
conf.lowRes = 0; % Downsampling the stack (1: yes, 0: no)
conf.dSamp = 1; % Downsampling amount
conf.numDil = 0; % Number of dilation
conf.numErd = 3; % Number of erosion
conf.in_cor_dm = 3; % inner core diameter
conf.out_cor_dm = 2.5; % Outer core diameter
conf.halfPortionChop = 0.3; % how much of nuclear slice will be excluded on the top (30%) and bottom (30%) for rim analyis
conf.paramsName = {'Timepoint' ,'Vol_Nuc1', 'SurfArea_Nuc1', 'Vol_Nuc2', 'SurfArea_Nuc2', 'AvgIntOuterCore1_ms', 'AvgIntInnerCore2_ms', 'AvgIntInnerCore3_ms', 'AvgIntOuterCore4_ms', ...
              'AvgIntNonCore12_ms', 'AvgIntNonCore34_ms', 'TotPixOuterCore1_ms', 'TotPixInnerCore2_ms', 'TotPixInnerCore3_ms', 'TotPixOuterCore4_ms', 'TotPixNonCore12_ms', ...
              'TotPixNonCore34_ms', 'AvgIntNuc1_ms', 'TotPixNuc1_ms', 'AvgIntNuc2_ms', 'TotPixNuc2_ms', 'AvgIntOuterCore1_ss', 'AvgIntInnerCore2_ss', 'AvgIntInnerCore3_ss', ...
              'AvgIntOuterCore4_ss', 'AvgIntNonCore12_ss', 'AvgIntNonCore34_ss', 'TotPixOuterCore1_ss',  'TotPixInnerCore2_ss', 'TotPixInnerCore3_ss', 'TotPixOuterCore4_ss', ...
              'TotPixNonCore12_ss', 'TotPixNonCore34_ss', 'AvgIntNuc1_ss', 'TotPixNuc1_IsoRim', 'AvgIntNuc2_ss', 'TotPixNuc2_IsoRim', 'Pix_ratio_lhi_vs_uhi'};
end