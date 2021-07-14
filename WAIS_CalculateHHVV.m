function HHVV = WAIS_CalculateHHVV(E,S,cfg) 
% Evaluates HHVV parameters. 
% 
% Input: 
% - E: Output of WAIS_ScatteringMatrix.
% - S: Output of WAIS_ReconstructAzimuth. 
% - cfg: File containing all configuration parameters.
%
% Output: 
% - Azimuth: Vector of downsampled azimuth [degrees]
% - Depth: Vector of downsampled depth [m]
% - C: HHVV coherence, range 0 (no coherence) to 1 (full coherence)
% - P: HHVV phase angle [rad]
% - Pu: HHVV phase angle, unwrapped [rad]
% - Pe: HHVV phase angle error (std) [rad]
% - dP: HHVV phase derivative through depth [rad m^-1]
% - dP_u: Upper bound of HHVV phase derative
% - dP_l: Lower bound of HHVV phase derative
%
% TJ Young
% 06 February 2020

%% Prep work

% Re-label commonly used variables
Azimuth = E.Azimuth;
Depth = E.Depth;
Shh = S.Shh; 
Svv = S.Svv; 

% Calculate resolutions
dx = mean(diff(Azimuth * pi/180)); % Angle resolution [rad]
dz = mean(diff(Depth)); % Depth resolution [m]

% Convert window size to bins
xWin = cfg.xWin * pi/180; % Azimuth window size (radians)
xWin = round(xWin/dx); % Convert azimuth window size to bins
zWin = round(cfg.zWin/dz); % Convert depth window size to bins

% Vectors to arrays 
Azimuths = repmat(Azimuth,length(Depth),1)';
Depths = repmat(Depth,length(Azimuth),1);

%% Calculate polarimetric coherence (C_hhvv)

clear Depth_hhvv Azimuth_hhvv C_hhvv
zWinIdx = (1:1:zWin) - zWin; % Initiate z-moving window
for ii = 1:floor(length(Depth)/zWin)
    xWinIdx = (1:1:cfg.xWin) - cfg.xWin; % Initiate x-moving window
    zWinIdx = zWinIdx + zWin; % Move z-window
    Depth_hhvv(ii) = mean(Depth(zWinIdx)); % Downsample depth vector
    for jj = 1:floor(length(Azimuth)/cfg.xWin)
        xWinIdx = xWinIdx + cfg.xWin; % Move x-window
        Azimuth_hhvv(jj) = mean(Azimuth(xWinIdx)); % Downsample azimuth vector
        %Numer = sum(Shh(xWinIdx,zWinIdx) .* conj(Svv(xWinIdx,zWinIdx)),'all');
        Numer = sum(conj(Shh(xWinIdx,zWinIdx)) .* Svv(xWinIdx,zWinIdx),'all');
        Denom1 = sqrt(sum(abs(Shh(xWinIdx,zWinIdx)).^2,'all'));
        Denom2 = sqrt(sum(abs(Svv(xWinIdx,zWinIdx)).^2,'all'));
        C_hhvv(jj,ii) = Numer ./ (Denom1.*Denom2); % Eq. 19 of Jordan et al. (2019)
        C_hhvv(jj,ii) = conj(C_hhvv(jj,ii)); % Obtain deramped phase, Jordan et al. (2020) Eq. 7
    end
end

P_hhvv = angle(C_hhvv); % hhvv phase difference, Jordan et al. (2019) Eq. 14
P_hhvv_u = unwrap(P_hhvv,[],2); % hhvv phase difference unwrapped

% Calculate phase error of C_hhvv
CramerRao = @(CHHVV,N) (1./abs(CHHVV)) .* sqrt((1-abs(CHHVV).^2)./(2.*N));
Pe_hhvv = CramerRao(C_hhvv,cfg.zWin);

% Calculate lower and upper bounds of C_hhvv
C_hhvv_upper = C_hhvv .* complex(cos(Pe_hhvv),sin(Pe_hhvv));
C_hhvv_lower = C_hhvv .* complex(cos(-Pe_hhvv),sin(-Pe_hhvv));

%% Calculate hhvv variables

dx_hhvv = mean(diff(Azimuth_hhvv)); % Re-calculate dx for hhvv downsampled resolution
dz_hhvv = mean(diff(Depth_hhvv)); % Re-calculate dz for hhvv downsampled resolution

% Calculate hhvv phase gradient (dP_hhvv), Jordan et al. (2019) Eq. 23
DP_hhvv = RI(C_hhvv,dz_hhvv,cfg);
DP_hhvv_u = RI(C_hhvv_upper,dz_hhvv,cfg);
DP_hhvv_l = RI(C_hhvv_lower,dz_hhvv,cfg);

%% Export variables to structure

HHVV = struct;
HHVV.Azimuth = Azimuth_hhvv; 
HHVV.Depth = Depth_hhvv;
HHVV.C = C_hhvv; 
HHVV.P = P_hhvv; 
HHVV.Pu = P_hhvv_u; 
HHVV.Pe = Pe_hhvv;
HHVV.dP = DP_hhvv;
HHVV.dP_u = DP_hhvv_u;
HHVV.dP_l = DP_hhvv_l; 

end

%% Calculate P_hhvv using real and imaginary components 
function RIout = RI(CHHVV,dZ,cfg)

%global cfg

Real = real(CHHVV);
Imag = imag(CHHVV);
%Real = pkConvol(Real,cfg.ftype,cfg.fsize,cfg.fsigma,cfg.fwindow);
%Imag = pkConvol(Imag,cfg.ftype,cfg.fsize,cfg.fsigma,cfg.fwindow);
Real = imgaussfilt(Real,cfg.fsigma,'FilterSize',cfg.fwindow);
Imag = imgaussfilt(Imag,cfg.fsigma,'FilterSize',cfg.fwindow);


for ii = 1:size(Real,1)
    dRdz(ii,:) = gradient(Real(ii,:),dZ);
    dIdz(ii,:) = gradient(Imag(ii,:),dZ);
end

RIout = (Real.*dIdz - Imag.*dRdz) ./ (Real.^2 + Imag.^2);

end