function S = WAIS_ReconstructAzimuth(E,cfg)
% Evaluates quad-pol reconstruction using E.
% 
% Input: 
% - E: Output of WAIS_ScatteringMatrix.
% - cfg: File containing all configuration parameters.
%
% Output: 
% - S: Reconstructed quad-pol arrays
% - PN: Normalised power of HH and VH
% - PD: Phase difference of HH and VH
% - AzimuthMid: Mid-point of Azimuth vector (as PD is length(Azimuth)-1)
%
% TJ Young (originally from Carlos Martin)
% 06 February 2020

%% Re-label commonly used variables

Ehh = E.Ehh;
Ehv = E.Ehv;
Evh = E.Evh;
Evv = E.Evv;
Azimuth = E.Azimuth;
Depth = E.Depth;

%% Evaluate Fujita experiment

% Pre-allocate matrices
Shh = nan(cfg.nShots,length(Ehh)); 
Shv = nan(cfg.nShots,length(Ehv)); 
Svh = nan(cfg.nShots,length(Evh)); 
Svv = nan(cfg.nShots,length(Evv)); 
PhaseDiffPar = nan(cfg.nShots-1,length(Ehh));
PhaseDiffPer = nan(cfg.nShots-1,length(Ehv));

% Run Fujita model (Eq. 4 of Young et al. 2020)
Disp(['Applying Fujita et al. (2006) model...'])
for i=1:cfg.nShots
    t = Azimuth(i) * (pi/180); % Degrees to radians
    Shh(i,:)=Ehh*cos(t)^2 - (Ehv+Evh)*cos(t)*sin(t) + Evv*sin(t)^2; % Originally ErPar
    Svv(i,:)=Evv*cos(t)^2 + (Ehv+Evh)*cos(t)*sin(t) + Ehh*sin(t)^2;
    Shv(i,:)=(Ehh-Evv)*cos(t)*sin(t) + Ehv*cos(t)^2 - Evh*sin(t)^2; % Originally ErPer
    Svh(i,:)=(Ehh-Evv)*cos(t)*sin(t) - Ehv*sin(t)^2 + Evh*cos(t)^2;
end

%% Calculate power anomaly and phase difference
% From Brisbourne et al. (2019) 
% Power anomaly from Eq. 3 of Matsuoka et al. (2012) 
% Phase difference from Eq. 5 of Young et al. (2020)

Disp(['Calculating power anomaly...'])

PowerPar=20*log10(abs(Shh));
PowerNormPar=detrend(PowerPar,'constant');
PowerNormPar=AverageDepth(PowerNormPar',Depth,cfg.dzPowerNorm)';

PowerPer=20*log10(abs(Svh));
PowerNormPer=detrend(PowerPer,'constant');
PowerNormPer=AverageDepth(PowerNormPer',Depth,cfg.dzPowerNorm)';

Disp(['Calculating phase difference...'])

AzimuthMid = (Azimuth(1:end-1)+Azimuth(2:end))/2;

for j=1:cfg.nShots-1 % Clockwise or counterclockwise 
     PhaseDiffPar(j,:)=angle(conj(Shh(j+1,:)).*Shh(j,:));
     PhaseDiffPer(j,:)=angle(conj(Svh(j+1,:)).*Svh(j,:));
end

xdiffphase=AverageDepth(cos(PhaseDiffPar'),Depth,cfg.dzdiffphase);
ydiffphase=AverageDepth(sin(PhaseDiffPar'),Depth,cfg.dzdiffphase);
PhaseDiffPar=atan2(ydiffphase,xdiffphase)';

xdiffphase=AverageDepth(cos(PhaseDiffPer'),Depth,cfg.dzdiffphase);
ydiffphase=AverageDepth(sin(PhaseDiffPer'),Depth,cfg.dzdiffphase);
PhaseDiffPer=atan2(ydiffphase,xdiffphase)';

%% Export variables to structure

S = struct;
S.Shh = Shh; 
S.Shv = Shv;
S.Svh = Svh;
S.Svv = Svv; 
S.PNPar = PowerNormPar;
S.PNPer = PowerNormPer; 
S.PDPar = PhaseDiffPar;
S.PDPer = PhaseDiffPer; 
S.AzimuthMid = AzimuthMid; 

