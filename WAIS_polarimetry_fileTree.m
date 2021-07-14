function [mdF,fileNames,I] = WAIS_polarimetry_fileTree(Point)
% Automating processing and plotting of polarimetry through depth of sites
% along a transect.
% 
% File 'WAIS_metadata.xlsx' can be obtained from Young and Dawson (2021)
% via the BAS Polar Data Centre.
% 
% Note: Make sure to change the inFolder variable to point to the location
%   of the script repository.
% 
% References:
% Young, T. J. and Dawson, E. J. (2021) Quad-polarimetric ApRES
%   measurements  along a 6 km-long transect at the WAIS Divide, December 
%   2019 (Version 1.0) [Data set], NERC EDS UK Polar Data Centre, 
%   doi:10.5285/BA1CAF7A-D4E0-4671-972A-E567A25CCD2C.
%
% TJ Young (originally from Carlos Martin)
% 06 February 2020

%% A. Set up script and file structure
%inFolder = ('/Volumes/GoogleDrive/My Drive/Manuscripts/2020Young_WAISD/BAS_PDC/Young_WAISDivide_ApRES_polarimetry/');
inFolder = ('H:/My Drive/Manuscripts/2020Young_WAISD/BAS_PDC/Young_WAISDivide_ApRES_polarimetry/');
addpath(genpath(inFolder));

%% B. Obtain file names from metadata

% Read in metadata and reformat columns
inFile = 'WAIS_metadata.xlsx';
opts = detectImportOptions(inFile);
opts.VariableTypes(5) = {'double'};
md = readtable(inFile,opts,'ReadVariableNames',false); 
md.Site = categorical(md.Site);
md.Orientation = categorical(md.Orientation);
md.Direction = categorical(md.Direction);
md.Mode = categorical(md.Mode); 
md.Date_ = datestr(md.Date_,'yyyy-mm-dd');
md.Time_ = datestr(md.Time_,'HHMMSS');

% Find index of rows matching point
mdF = md(Point == md.Site,:);

% Reconstruct filenames from metadata
for ii = 1:4
    if mdF.Mode(ii) == 'A'
        fFN = strcat('WAISD_',Point,'_',cellstr(mdF.Orientation(ii)),...
            '_',mdF.Date_(ii,:),'_',mdF.Time_(ii,:),'.dat');
    elseif mdF.Mode(ii) == 'U'
        fFN = strcat('WAISD_',Point,'_',cellstr(mdF.Orientation(ii)),...
            '_DATA',mdF.Date_(ii,:),'-',mdF.Time_(ii,1:4),'.DAT');
    end
    fileNames(ii,1) = fFN;
end

% Index of orientations
I.vv = find(mdF.Orientation == 'VV');
I.vh = find(mdF.Orientation == 'VH');
I.hh = find(mdF.Orientation == 'HH');
I.hv = find(mdF.Orientation == 'HV');