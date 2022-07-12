%%
% This package of codes quantifies the distributions of Nups insdie the core
% and non-core regions on the nuclear surface. This script is main funtion
% that gets configurations (running options and parameters) and call the 
% 'fn_quantify_nups_kinetics' function to batch process all the tif files
% inside the input directory.
% Author: M. Julius Hossain, EMBL, Heidelberg
%%

% Add directories containing files supporting input output
addpath 'IO'
addpath(fullfile('IO', 'bioformats'));

% Reads the configuration file
conf = fn_configurations();
% Channel index of chromatin
chr_chan_idx = conf.chr_chan_idx;
% Channel index of nup
prot_chan_idx = conf.prot_chan_idx;

% Root directory containing the data
fp = conf.in_dir;

% Get the list of filenames inside the data directory
fn = dir(fullfile(fp, '*.tif'));

% Process individual image stack one by one
for file_idx = 1:length(fn)
    try
        fn_quantify_nups_kinetics(fp, fn(file_idx).name, chr_chan_idx, prot_chan_idx, file_idx, conf);
    catch
        % Display the filename if it could not be processed and move on to
        % the next tif file.
        disp(['Could not process: ' fn(file_idx).name]);
    end
end
