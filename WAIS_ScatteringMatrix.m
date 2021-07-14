function E = WAIS_ScatteringMatrix(fileNames,cfg)
% Load quad-pol measurements & conduct initial pre-processing.
% 
% Input: 
% - fileNames: 4x1 cell array of file names. Default settings assume order [HH;HV;VH;VV].
% - cfg: File containing all configuration parameters.
%
% Output: 
% - E: Complex vector of SpecCor of each file
% - Azimuth: Vector of all azimuth records [degrees]
% - Depth: Vector of all depth records [m]
%
% TJ Young
% 06 February 2020

%% Load data

% Determine index of orientations
if length(cfg.fileOrder) == 4 % If order is determined in config
    I.hh = cfg.fileOrder(1); 
    I.hv = cfg.fileOrder(2);
    I.vh = cfg.fileOrder(3);
    I.vv = cfg.fileOrder(4);
else 
    try % Determine order if it is embedded in file names
        I.hh = find(contains(fileNames,'_HH_'),1);
        I.hv = find(contains(fileNames,'_HV_'),1);
        I.vh = find(contains(fileNames,'_VH_'),1);
        I.vv = find(contains(fileNames,'_VV_'),1);
    catch % Else, use default order
        I.hh = 1; I.hv = 2; I.vh = 3; I.vv = 4; 
    end
end

Disp('Loading quad-polarimetric files...')
cfg.ei = mean([cfg.eiPar,cfg.eiPer]);
[~,fhh] = evalc('fmcw_load(fileNames{I.hh},1,cfg.ei);');
[~,fhv] = evalc('fmcw_load(fileNames{I.hv},1,cfg.ei);');
[~,fvh] = evalc('fmcw_load(fileNames{I.vh},1,cfg.ei);');
[~,fvv] = evalc('fmcw_load(fileNames{I.vv},1,cfg.ei);');

%% Clean data

% Manually select chirps
if cfg.doManualChirpSelect
    disp('Using chirp subset defined in config')
    fvv = fmcw_burst_subset(fvv,cfg.chirplist);
    fvh = fmcw_burst_subset(fvh,cfg.chirplist);
    fhh = fmcw_burst_subset(fhh,cfg.chirplist);
    fhv = fmcw_burst_subset(fhv,cfg.chirplist);
end

% Frequency range
if cfg.fRange(1) > 2e8 || cfg.fRange(2) < 4e8
    % Cull f.specCor range
    disp(['Cropping frequency range to: ' num2str(cfg.fRange(1)/1e6) '-' num2str(cfg.fRange(2)/1e6) ' MHz'])
    disp(' ')
    fvv = fmcw_cull_freq(fvv,cfg.fRange);
    fvh = fmcw_cull_freq(fvh,cfg.fRange);
    fhh = fmcw_cull_freq(fhh,cfg.fRange);
    fhv = fmcw_cull_freq(fhv,cfg.fRange);
end

% Clean shots
if cfg.cullContaminated
    % Cull gross contaminated chirps
    [fvv,nBad] = fmcw_cull_bad(fvv,cfg.noisePowerLimit,0);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' fileNames{I.vv}])
    [fvh,nBad] = fmcw_cull_bad(fvh,cfg.noisePowerLimit,0);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' fileNames{I.vh}])
    [fhh,nBad] = fmcw_cull_bad(fhh,cfg.noisePowerLimit,0);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' fileNames{I.hh}])
    [fhv,nBad] = fmcw_cull_bad(fhv,cfg.noisePowerLimit,0);
    Disp(['Removed ' int2str(nBad) ' contaminated chirps from ' fileNames{I.hv}])
end
    
% Cull noisy chirps
if cfg.cullNoisy
    Disp(['Culling noisiest ' int2str(cfg.nNoisest) ' chirps...'])
    fvv = fmcw_cull_noisey(fvv,cfg.nNoisest,0);
    fvh = fmcw_cull_noisey(fvh,cfg.nNoisest,0);
    fhh = fmcw_cull_noisey(fhh,cfg.nNoisest,0);
    fhv = fmcw_cull_noisey(fhv,cfg.nNoisest,0);
end

%% Average and migrate bursts (if necessary)

% Average all chirps to one sub-burst
fvv = fmcw_burst_mean(fvv);
fvh = fmcw_burst_mean(fvh);
fhh = fmcw_burst_mean(fhh);
fhv = fmcw_burst_mean(fhv);

% Migrate vif voltage if different attenuations have been used
fN = fileNames{1};
if contains(fN,{'WAISD'})
if ~contains(fN(7),{'A','H','I','J'})
    Disp(['Scaling VH | HV power...'])
    fvh.vif0 = fvh.vif; % Original vif of fvhm
    fhv.vif0 = fhv.vif; % Original vif of fhvm
    [fvh.vif,fhv.vif] = scaleVifByAtt(fvv.vif,fvv.Attenuator_1,...
        fvh.vif0,fvh.Attenuator_1,fhh.vif,fhh.Attenuator_1,...
        fhv.vif0,fhv.Attenuator_1);
end
end

%% Range process bursts

Disp('Range processing bursts...')
[fvv.rangeCoarse,fvv.rangeFine,fvv.specCor,fvv.specRaw] = fmcw_range(fvv,cfg.p,cfg.rangeLim(2),cfg.winFun);
[fvh.rangeCoarse,fvh.rangeFine,fvh.specCor,fvh.specRaw] = fmcw_range(fvh,cfg.p,cfg.rangeLim(2),cfg.winFun);
[fhh.rangeCoarse,fhh.rangeFine,fhh.specCor,fhh.specRaw] = fmcw_range(fhh,cfg.p,cfg.rangeLim(2),cfg.winFun);
[fhv.rangeCoarse,fhv.rangeFine,fhv.specCor,fhv.specRaw] = fmcw_range(fhv,cfg.p,cfg.rangeLim(2),cfg.winFun);
zi=fhh.rangeCoarse; % Coarse range
Azimuth = 0:cfg.dx:180;

%% Set antenna orientation

eval(['Evv = ',cfg.antOrient(4),'fvv.specCor;'])
eval(['Evh = ',cfg.antOrient(3),'fvh.specCor;'])
eval(['Ehh = ',cfg.antOrient(1),'fhh.specCor;'])
eval(['Ehv = ',cfg.antOrient(2),'fhv.specCor;'])

%% True depth correction
% From Carlos Martin

C.rhoi = 910; % Density of ice 
C.rhos = 400; 
C.Lrho = 40;
C.nI = 1.68;
C.c = 299792458; %speed of light

ShouldBeZero = @(d,dI,L,Rhosp) -dI+d+L*(C.nI-1)/C.nI*(1-Rhosp)*(exp(-d/L)-1);
TrueDepthfun = @(Rhosp,L,dI) bisection(@(d) ShouldBeZero(d,dI,L,Rhosp),0,max(dI)+20);
Depth = TrueDepthfun(C.rhos/C.rhoi,C.Lrho,zi);

%% Export variables to structure

E = struct; % Structure for scattering matrix elements
E.Ehh = Ehh; 
E.Ehv = Ehv;
E.Evv = Evv;
E.Evh = Evh; 
E.Azimuth = Azimuth; % Azimuthal resolution
E.Depth = Depth; % Ice depth

function [vifVHT,vifHVT] = scaleVifByAtt(vifVV,attVV,vifVH,attVH,vifHH,attHH,vifHV,attHV)

% Find mean and standard deviation of HV/VH shots that need migration
VV.mean = mean(vifVV);
VH.mean = mean(vifVH);
HH.mean = mean(vifHH);
HV.mean = mean(vifHV);
VV.preStd = std(vifVV); 
VH.preStd = std(vifVH); 
HH.preStd = std(vifHH); 
HV.preStd = std(vifHV); 

% Load parameters of linear regression (derived from script scaleViffByAttPre.m)
m.hi  = -0.013239142904390;
m.mid = -0.010636859641062;
m.lo  = -0.008034576377735;
b.hi  =  0.400527445258937;
b.mid =  0.317261567590820;
b.lo  =  0.233995689922704;
        
% Linear regression
y = @(m,x,b) m*x+b;

% Migrate subburst to appropriate power StD
preMean = y(m.mid,attHV,b.mid);
preStd = preMean - y(m.lo,attHV,b.lo);
postMean = y(m.mid,attHH,b.mid);
postStd = postMean - y(m.lo,attHH,b.lo);
postHi = postMean + postStd;
postLo = postMean - postStd;

% Find percentile of HV/VH shots taht need migration
HV.preStdPrct = normcdf(HV.preStd,preMean,preStd);
VH.preStdPrct = normcdf(VH.preStd,preMean,preStd);
% Find corresponding STD of post-migrated shot 
HV.postStd = norminv(HV.preStdPrct,postMean,postStd);
VH.postStd = norminv(VH.preStdPrct,postMean,postStd);
% Move to mean-std if postStd is too low 
if HV.postStd > postHi, HV.postStd = postHi; end
if VH.postStd > postLo, VH.postStd = postHi; end
if HV.postStd < postLo, HV.postStd = postLo; end
if VH.postStd < postLo, VH.postStd = postLo; end

% Migrate HV/VH shots 
vifHVT = (vifHV - HV.mean)/HV.preStd;
vifHVT = vifHVT .* HV.postStd;
vifHVT = vifHVT + HV.mean; 
vifVHT = (vifVH - VH.mean)/VH.preStd;
vifVHT = vifVHT .* VH.postStd;
vifVHT = vifVHT + VH.mean; 