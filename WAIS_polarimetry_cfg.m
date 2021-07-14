% Config parameters for all of TIME Polarimetry files
% 
% TJ Young
% 03 March 2020

global cfg

% GENERAL: Display
cfg.verbose = 1; % Display intermediary status reel
cfg.saveFigs = 1; % Save figures
cfg.figType = '.pdf'; % Figure output file type ['.pdf' '.png' '.fig']

% GENERAL: Setup and dimensions
cfg.fileOrder = 0; % Determine file order, either 0 (use default) or a four element vector (order [HH HV VH VV])
cfg.AzimuthRange = [-90 90]; % Range of Azimuths, degrees (either 0 to 180 or -90 to 90)
cfg.nShots=180+1; % Number of orientations, +1 to account for phase wrapping 
cfg.dx = (range(cfg.AzimuthRange)+1)/(cfg.nShots);
cfg.antOrient = '++--'; % Antenna orientation [order HH HV VH VV]

% Chirp cleaning
cfg.doManualChirpSelect = 0;
cfg.fRange = [2e8 4e8]; % Tx frequency range to use - full range
cfg.chirplist = [1:100];
cfg.cullContaminated = 0;
cfg.cullNoisy = 1;
cfg.noisePowerLimit = 0.005; % allowable percent difference in power from quadratic fit to chirpnum-power trend
cfg.nNoisest = 20;

% Range processing parameters
cfg.p = 2; % Pad factor
cfg.bulkAlignRange=[100 200]; % Bulk align range
cfg.rangeLim = [0 1500]; % max range to plot
cfg.maxRange = 3000;
cfg.winFun = @blackman; % @blackman @blackmanharris

% Power anomaly and phase differencing parameters
cfg.dzPower=10; % 
cfg.dzPowerNorm=10;
cfg.dzdiffphase=30;

% PROCESS: Switches
cfg.firncorrection = 1; % Calculate anisotropy parameters for firn correction
cfg.doFilter = 1; % Filter results by noise floor

% hhvv parameters
cfg.calculateE2E1 = 1; % Calculates E2-E1
cfg.xWin = 1; % Azimuth window size [degrees]
cfg.zWin = 15; % Depth window size for HHVV [m]
cfg.ftype = 'gaussian';
cfg.fsigma = 2.2; % Gaussian filter sigma
cfg.fwindow = [1 5]; % Gaussian filter window [bins]

% Fast axes parameters
cfg.assignE2E1 = 'user'; % Fast axes assignment method: 'user' 'auto'
cfg.FA0 = 90; % Fast axes (E1) 'user' assignment (scalar, or nx2 array of [azimuths depths])
cfg.doSmooth = 1; % Apply smoothing to fast axes
cfg.zWinA = 100; % Depth window size for fast axes smoothing [m]
cfg.prom = 2; % Fast Axes peak prominence

% Constants
cfg.fc = 3e8; % Centre frequency [Hz]
cfg.eiPar = 3.169; % Parallel component of real part of dielectric permittivity
cfg.eiPer = 3.134; % Perpendicular component of real part of dielectric permittivity
cfg.dei = 0.034; % Dielectric anisotropy
cfg.c = 299792458; % Speed of light [m a^-1]