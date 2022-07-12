# Source code for Nups intensity quantification for Otsuka et al. Nature_2022
This source code extracts distributions of nups intensity on nulear surface in core and non-core regions. This was partly developed for the manuscript "A quantitative map of nuclear pore assembly reveals two distinct mechanisms" by Otsuka et al. Nature 2022 which is under publication. A biorxiv version of the manuscript is available at doi: https://doi.org/10.1101/2021.05.17.444137
## Requirements
The code has been tested on the following versions. Previous versions of MATLAB may also work:  
* MATLAB 2022a  
* Bioformat for MATLAB (http://downloads.openmicroscopy.org/bio-formats/5.3.4/)
# How to run
* Set the configurations (running options and parameters) in 'fn_configurations.m' file
* Run 'fn_quantify_nups_kinetics_main'
* The configurations parameters, segmentation visualizations and extracted parameters are saved in the output folder
