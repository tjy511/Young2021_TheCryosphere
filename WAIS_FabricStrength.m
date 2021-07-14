function A = WAIS_FabricStrength(E,HHVV,X,cfg)
% Calculates orientation of E1 and E2 eigenvectors.
% 
% Input: 
% - E: Output of WAIS_ScatteringMatrix.
% - HHVV: Output of WAIS_CalculateHHVV.
% - X: Output of WAIS_FabricOrientation. 
% - cfg: File containing all configuration parameters.
%
% Output: 
% - a: Azimuthal fabric strength
% - a_se: Std of a
% - E1: Downsampled E1 to HHVV resolution
% - E2: Downsampled E2 to HHVV resolution
%
% TJ Young
% 06 February 2020

%% Re-label commonly used variables
Azimuth = E.Azimuth;
Depth = E.Depth;
E1 = X.E1;
E2 = X.E2;

%% Downsample fast axes to HHVV resolution

% Calculate bins to downsample 
dz = mean(diff(Depth)); % Depth resolution [m]
zWin = round(cfg.zWin/dz); 

zWinIdx = (1:1:zWin) - zWin; % Initiate z-moving window
for ii = 1:floor(length(Depth)/zWin)
    zWinIdx = zWinIdx + zWin; % Move z-window
    Depth_hhvv(ii) = mean(Depth(zWinIdx)); % Downsample depth vector
    
    % Indices within z-window
    ci0 = E1(zWinIdx); 
    ci90 = E2(zWinIdx); 
    
    % Mean values within z-window
    E1_hhvv(ii) = meanAngle180(ci0); 
    E2_hhvv(ii) = meanAngle180(ci90); 
    
    % Find fast axes index within HHVV azimuth
    [~,E1_hhvvi(ii)] = closest(HHVV.Azimuth,E1_hhvv(ii));
    [~,E2_hhvvi(ii)] = closest(HHVV.Azimuth,E2_hhvv(ii));
    
    % Extract dP_hhvv/dz values from fast axes indices
    dP1(ii) = HHVV.dP(E1_hhvvi(ii),ii);
    dP2(ii) = HHVV.dP(E2_hhvvi(ii),ii);
    dP1_u(ii) = HHVV.dP_u(E1_hhvvi(ii),ii);
    dP1_l(ii) = HHVV.dP_l(E1_hhvvi(ii),ii);
    dP2_u(ii) = HHVV.dP_u(E2_hhvvi(ii),ii);
    dP2_l(ii) = HHVV.dP_l(E2_hhvvi(ii),ii);
end


%% Calculate E2E1

% Jordan et al. (2019) Eq. 22
CalculateE2E1 = @(P,ei,dei,c,fc,firnc) P .* (2.*sqrt(ei)./(dei.*firnc)) .* c./(4.*pi.*fc);

% Implement firn correction
if cfg.firncorrection
    firnC = WAIS_firn_correction(HHVV.Depth,cfg.eiPer);
else
    firnC = 1;
end

% Calculate azimuthal fabric asymmetry
ei = mean([cfg.eiPar,cfg.eiPer]);
if cfg.calculateE2E1
    a1 = CalculateE2E1(dP1,ei,cfg.dei,cfg.c,cfg.fc,firnC);
    a2 = CalculateE2E1(dP2,ei,cfg.dei,cfg.c,cfg.fc,firnC);
    a1_u = CalculateE2E1(dP1_u,ei,cfg.dei,cfg.c,cfg.fc,firnC);
    a1_l = CalculateE2E1(dP1_l,ei,cfg.dei,cfg.c,cfg.fc,firnC);
    a2_u = CalculateE2E1(dP2_u,ei,cfg.dei,cfg.c,cfg.fc,firnC);
    a2_l = CalculateE2E1(dP2_l,ei,cfg.dei,cfg.c,cfg.fc,firnC);
    
    A1_se = abs(a1_u-a1_l); 
    A2_se = abs(a2_u-a2_l);
    
    a = nanmean([abs(a1);abs(a2)],1);
    a_se = sqrt(A1_se.^2 + A2_se.^2);
end

%% Export variables to structure

A = struct;
A.a = a; 
A.a_se = a_se; 
A.E1 = E1_hhvv;
A.E2 = E2_hhvv; 

end 

%% Mean angle accounting for phase wrapping
function mu = meanAngle180(X,dim)

% Set default operating dimensions
if nargin < 2
    dim = 1; 
end

% Override default dimensions if X is a vector
if size(X,1) == 1
    dim = 2; 
end

XR = X * (pi/180*2); % Convert to radians and scale range to [0 360]
Xe = exp(1i*XR); % Euler's Formula
Xm = nansum(Xe,dim); % Average across imaginary plane (nanmean works too, but nansum is mathematically more rigorous)
mu = atan2(imag(Xm),real(Xm)) * (180/pi/2); % Convert back to degrees (calculate arg), rescale to [-90 +90]
mu = mod(mu,180); % Circ shift to [0 180]

end