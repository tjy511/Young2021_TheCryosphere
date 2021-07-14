function X = WAIS_FabricOrientation(E,S,HHVV,cfg)
% Calculates orientation of E1 and E2 eigenvectors.
% 
% Input: 
% - E: Output of WAIS_ScatteringMatrix.
% - S: Output of WAIS_ReconstructAzimuth. 
% - HHVV: Output of WAIS_CalculateHHVV. 
% - cfg: File containing all configuration parameters.
%
% Output: 
% - E1: Relative orientation of E1 eigenvector [degrees]
% - E2: Relative orientation of E2 eigenvector [degrees]
%
% TJ Young
% 06 February 2020

%% Prep work

% Re-label commonly used variables
Shh = S.Shh; 
Svh = S.Svh;  
Azimuth = E.Azimuth;
Depth = E.Depth;
Azimuths = repmat(Azimuth',[1 size(Svh,2)]);
Depths = repmat(Depth,[size(Svh,1) 1]);


%% Find symmetry axes
% Note that this section does not yet distinguish between E1 and E2

% Duplicate matrices to full 360 degrees (no wrapping)
Azimuth359 = [Azimuth(1:end-1) Azimuth(1:end-1)+180];
Azimuths359 = repmat(Azimuth359',[1 size(Svh,2)]);
Svh359 = [Svh(1:end-1,:);Svh(1:end-1,:)]; 

% Re-apply power anomaly algorithm to Svh
SvhNorm359 = 20*log10(abs(Svh359)); 
SvhNorm359=detrend(SvhNorm359,'constant');
SvhNorm359=AverageDepth(SvhNorm359',Depth,cfg.dzPowerNorm)';
SvhNorm359Smooth = smoothdata(SvhNorm359,1,'gaussian'); % Smooth data to remove jenks

% Find indices of local minima through depth
TF = islocalmin(SvhNorm359Smooth,1,'FlatSelection','center',...
        'MinSeparation',(round(size(SvhNorm359Smooth,1)-1)/4) * 0.75,...
        'MinProminence',cfg.prom);
TFi = false(size(Svh));
for ii = 1:length(Depth)
    TFxs = Azimuth359(TF(:,ii)'); % Find degree value of indices
    TFxs = mod(TFxs,180); % Wrap values over 180 degrees
    TFxs = round(TFxs); % Remove rounding errors
    TFxs = unique(TFxs); % Remove duplicates
    TFxs = TFxs + Azimuth(1); % Shift minima from [0 179] to azimuth domain
    if length(TFxs) == 2 % Only work with robust picks (i.e. there are only two identified minima)
        TFi(:,ii) = ismember(Azimuth,TFxs); % Obtain index (back in 181 degrees)
        SA(:,ii) = TFxs; % Angles of minima 
    else
        SA(:,ii) = [NaN;NaN];
    end
end

%% Assign X2 and X1 eigenvectors

switch cfg.assignE2E1
    case 'user' % User-input vector of symmetry axis
        if (size(cfg.FA0,1) == 1 && size(cfg.FA0,2) == 1)
            SA0 = repmat(cfg.FA0,[1 length(Depth)]);
            SA90 = repmat(cfg.FA0+90,[1 length(Depth)]);
        elseif size(cfg.FA0,2) == 2
            % Fill this out later...
        end
        SAmod = mod(SA,180); % Wrap SA to [0 180]; 
        
        % Assign vectors
        for ii = 1:length(Depths)
            
            % Indices that fall into remit of SA0
            inds0 = SA0(ii)-44:1:SA0(ii)+44; 
            inds90 = SA90(ii)-44:1:SA90(ii)+44;
            inds0 = mod(inds0,180); inds90 = mod(inds90,180); % Account for azimuth wrapping
            
            % Indices of SA that fall within SA0 remit
            SAmod0 = double(ismember(SAmod(:,ii),inds0)); 
            SAmod90 = double(ismember(SAmod(:,ii),inds90)); 
            SAmod0(SAmod0==0) = NaN; SAmod90(SAmod90==0) = NaN; 
            
            % Extract and collapse to value
            E1(ii) = nanmean(SAmod(:,ii) .* SAmod0,1); % a2 --> - dP_hhvv/dz
            E2(ii) = nanmean(SAmod(:,ii) .* SAmod90,1); % a1 --> + dP_hhvv/dz 
        end
        
    case 'auto' % Automated assignment based on polarity of dPhhvv/dz
        
        % Interpolate dP_hhvv/dz to original resolution
        HHVVx = repmat(HHVV.Azimuth',[1 size(HHVV.dP,2)]);
        HHVVz = repmat(HHVV.Depth,[size(HHVV.dP,1) 1]);
        dP = interp2(HHVVz,HHVVx,HHVV.dP,...
            Depths,Azimuths,'spline',NaN); 
        
        % Assign eigenvectors to symmetry axes
        for ii = 1:length(Depth)
            try
            dPi = dP(TFi(:,ii),ii); % Find values of dP at symmetry axes indices
            xi = Azimuth(TFi(:,ii));
            if dPi(1) < 0
                E1(ii) = xi(1); % x1 corresponds to - dP_hhvv/dz
                E2(ii) = xi(2); % x2 corresponds to + dP_hhvv/dz
            end
            if dPi(1) > 0
                E1(ii) = xi(2); % x1 corresponds to - dP_hhvv/dz
                E2(ii) = xi(1); % x2 corresponds to + dP_hhvv/dz
            end
            catch
                E1(ii) = NaN;
                E2(ii) = NaN;
            end
        end

end

%% Export variables to structure

X = struct;
X.E1 = E1;
X.E2 = E2;

