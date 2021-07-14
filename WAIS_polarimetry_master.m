% Automating processing and plotting of polarimetry through depth of sites
% along a transect. This compilation of scripts can be used to produce the
% data presented in Figures 3, 4, and 6 in Young et al. (2020). 
%
% The complete WAIS Divide polarimetric ApRES dataset can be downloaded
% from the UK Polar Data Centre: 
% https://doi.org/10.5285/BA1CAF7A-D4E0-4671-972A-E567A25CCD2C
% 
% Note 1: These scripts are dependent on the fmcw package from Stewart
%   (2018), which the user will need to obtain separately from Dr. Craig L. 
%   Stewart (NIWA, New Zealand). 
% Note 2: Users will need to edit various paths pointing to data locations
%   in WAIS_polarimetry_fileTree.m and WAIS_firn_correction.m. 
% 
% References
% Stewart, C. L. (2018) Ice-ocean interactions under the north-western 
%   Ross Ice Shelf, Antarctica (PhD thesis), University of Cambridge.
%   doi:10.17863/CAM.21483. 
% Young, T. J. and Dawson, E. J. (2021) Quad-polarimetric ApRES
%   measurements along a 6 km-long transect at the WAIS Divide, December 
%   2019 (Version 1.0) [Data set]. NERC EDS UK Polar Data Centre. 
%   doi:10.5285/BA1CAF7A-D4E0-4671-972A-E567A25CCD2C
% Young, T. J. et al. (2020) Rapid and accurate polarimetric radar 
%   measurements of ice crystal fabric orientation at the Western Antarctic
%   Ice Sheet (WAIS) Divide ice core site. The Cryosphere Discussions,
%   1-27. doi:10.5194/tc-2020-264. 
% 
% TJ Young
% 06 February 2020

%% Set up script and file structure

% Obtain file tree
Site = 'WAISD';
Point = 'I';

% Load file tree
% Will need to edit the inFolder variable to point to dataset location
[~,fileNames,~] = WAIS_polarimetry_fileTree(Point);

% Load config parameters
WAIS_polarimetry_cfg;

%% Run processing algorithms

% 1. Load quad-pol measurement & reconstruct azimuthal response
E = WAIS_ScatteringMatrix(fileNames,cfg); 

% 2. Reconstruct azimuthal orientations
S = WAIS_ReconstructAzimuth(E,cfg); 

% 3. Calculate HHVV parameters
HHVV = WAIS_CalculateHHVV(E,S,cfg); 

% 4. Determine eigenvector orientations
X = WAIS_FabricOrientation(E,S,HHVV,cfg);

% 5. Calculate azimuthal fabric asymmetry 
A = WAIS_FabricStrength(E,HHVV,X,cfg); 

% 6. Plot output results
F = WAIS_plotFabricResults(E,S,HHVV,X,A,cfg);

%% 3. Export figures
% This section requires the export_fig functionality: 
% https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig

if cfg.saveFigs
    outFolder = ('/Users/tjy511/Downloads/');
    outFile = [Site,'_',Point,'_Results'];
    gFormat = cfg.figType;
    set(F,'Color','w')
    export_fig(figHandle(ii),[outFolder,outFile,gFormat],'-q101')
end