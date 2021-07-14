
% Ice fabric model that produces predicted pRES output from given fabric
% types and anisotropy/birefringence parameters. 

% Note 1: Users will need to change the folderOut and fileOut variables to
%   match the intended figure output path. Figure saving is switched by the
%   saveData variable. 
% Note 2: Users will need to change the path to the location of WAIS Divide
%   ice core fabric measurements. This dataset can be downloaded from 
%   https://nsidc.org/data/NSIDC-0605/versions/1.

% References
% Fitzpatrick, J. J. et al. (2014) Physical properties of the WAIS divide
%   ice core. Journal of Glaciology 60(224):1140-1154.
%   doi:10.3189/2014JoG14J100
% Fujita, S. et al. (2006) Radio-wave depolarization and scattering within
%   ice sheets: A matrix-based model to link radar and ice-core 
%   measurements and its application. Journal of Glaciology 52(178):407-424.
%   doi:10.3189/172756506781828548
% Jordan, T. M. et al. (2019) A polarimetric coherence method to determine
%   ice crystal orientation fabric from radar sounding: Application to the
%   NEEM ice core region. IEEE Transactions on Geoscience and Remote 
%   Sensing 57(11):8641-8657. doi:10.1109/TGRS.2019.2921980
% Young, T. J. et al. (2020) Rapid and accurate polarimetric radar 
%   measurements of ice crystal fabric orientation at the Western Antarctic
%   Ice Sheet (WAIS) Divide ice core site. The Cryosphere Discussions,
%   1-27. doi:10.5194/tc-2020-264. 

% TJ Young (original script from Carlos Martin)
% V1.0 - 20 April 2020

%% User-defined parameters
% Parameters that matter:
% Fabric type
% A: Anisotropic scattering
% H: Ice thickness
% SA: Symmetry of axis
% X: Azimuthal resolution
% Z: Depth resolution
% ei: Real part of dielectric permittivity (parallel or perpendicular)

% Model parameters
fabricType = 'WAIS Divide Ice Core';
A = 5; % Anisotropic scattering [dB] (10 in Fujita et al. 2006)
H = 1800; % Ice thickness
SA = 89; % Axis of symmetry (azimuth of power anomaly minima)
X = 91; % Azimuthal resolution (number of shots)
Z = 1000; % Depth resolution (number of numerical layers)
dzPowerNorm = 50; % Moving average depth in metres

cLim = 5; % C-axis limits
saveData = 0;

% Wave transmission
eiPar = 3.169; % Parallel component of real part of dielectric permittivity
eiPer = 3.134; % Perpendicular component of real part of dielectric permittivity

%% Set up model structure

switch fabricType
    case 'Isotropic'
        E1g = 1/3;
        E2g = 1/3;
        E3g = 1/3;
    case 'Vertical Partial Girdle'
        E1g = 0;
        E2g = 0.2;
        E3g = 1-E2g;
    case 'Vertical Pole'
        E1g = 0;
        E2g = 0;
        E3g = 1;
    case 'Vertical Girdle'
        E1g = 0;
        E2g = 1/2;
        E3g = 1/2;
    case 'Horizontal Pole'
        E1g = 0;
        E2g = 1;
        E3g = 0;
    case 'Horizontal Girdle'
        E1g = 1/2;
        E2g = 1/2;
        E3g = 0;
    case 'WAIS Divide Ice Core'
        % Import WAIS Divide Ice Core Fabric Eigenvalues
        % Note that ice core eigenvalues are S1 > S2 > S3
        %inFolder = '/Volumes/GoogleDrive/My Drive/ITGC-TIME/Datasets/WAIS Divide Ice Core Fabric/';
        inFolder = 'H:/My Drive/ITGC-TIME/Datasets/WAIS Divide Ice Core Fabric/';
        inFile = 'Database of Analyzed Sections.xls';
        opts = detectImportOptions([inFolder,inFile]);
        WF = readtable([inFolder,inFile],'Sheet','Sheet1');
        WF.Properties.VariableNames{6} = 'S2';
        
        % Overwrite model parameters
        H = WF.Section';
        
        % Define eigenvalues
        E1g = WF.S3; E1g = [1/3;E1g]';
        E2g = WF.S2; E2g = [1/3;E2g]';
        E3g = WF.S1; E3g = [1/3;E3g]';
        
        % Re-define anisotropic scattering
        Ax = [repmat(1,1,150),1:3/550:4,4:2/700:6,...
             repmat(6,1,2001)];
        A = downsample(Ax,round(length(Ax)/length(E1g)));
        
end

%% Set model parameters and constants

% Constants
temp = -20; % Temperature, in Celsius
e_0 = 8.8541878128e-12; % Electric permittivity in vacuum [F m^-1]
mu_0 = 4e-7*pi; % Magnetic permeability in vacuum [H m^-1]
c = 299792458; % Speed of light [m a^-1]

% Radar parameters
fc = 3e8; % Central frequency [Hz]
omega = 2*pi*fc; % Angular frequency

% Define geometry
zg = [0,H];

% Resize eigenvalues
if size(E1g) ~= size(zg)
    E1g = repmat(E1g,size(zg));
    E2g = repmat(E2g,size(zg));
    E3g = repmat(E3g,size(zg));
end

% Fabric alignment
phase0g=repmat(SA,size(zg)); % Fast axes directions

% Wave transmission
ei = mean([eiPar eiPer]);
dei = eiPar-eiPer; % Dielectric anisotropy

% Scattering parameters
Sxg = repmat(1,size(zg));
Syg = Sxg .* A;

% Prepare model input

disp(['Modelling fabric type: ', fabricType])

% Interpolate to numerical domain
z=linspace(0,H(end),Z);
phase0i=interp1(zg,phase0g,z)*pi/180;
Sxi=interp1(zg,Sxg,z);
Syi=interp1(zg,Syg,z);

% Firn correction parameter
if strcmp(fabricType,'WAIS Divide Ice Core')
    firnCi = WAIS_firn_correction(z,eiPer);
else
    firnCi = ones(size(z));
end

% Transform fabric eigenvalues to dielectric tensor
% Eq. A3 of Fujita et al. (2006)
% Here eigenvalues are using convention (E3 > E2 > E1)
exg=eiPer*ones(size(zg))+E1g*dei;
eyg=eiPer*ones(size(zg))+E2g*dei;
exi=interp1(zg,exg,z);
eyi=interp1(zg,eyg,z);% .* firnCi;

% Calculate average on the layer
ex=(exi(1:Z-1)+exi(2:Z))/2;
ey=(eyi(1:Z-1)+eyi(2:Z))/2;
phase0=(phase0i(1:Z-1)+phase0i(2:Z))/2;
phase0deg = rad2deg(phase0);
Sx=(Sxi(1:Z-1)+Sxi(2:Z))/2;
Sy=(Syi(1:Z-1)+Syi(2:Z))/2;

% Components of T-Matrix
k0 = 2*pi/(c/fc); % Wave number in a vacuum
kx=2*pi*fc/c*sqrt(ex); % Propagation speed along x
ky=2*pi*fc/c*sqrt(ey); % Propagation speed along y
dz=z(2:end)-z(1:end-1);
T_ix = exp(-1i.*k0.*dz + 1i.*kx.*dz); % Eq. 6a of Fujita et al. (2006)
T_iy = exp(-1i.*k0.*dz + 1i.*ky.*dz); % Eq. 6a of Fujita et al. (2006)

ne=Z-1; % Number of elements
R=zeros(2,2,ne); %Rotation Matrix
T=zeros(2,2,ne); %transmission matrix

% Set domain
x = 0:pi/(X-1):pi;
xdeg = x * 180/pi;
dx = x(2:end)-x(1:end-1); 

%% Run Fujita et al. (2006) model

for j=1:X
    
    % Establish matrices
    for i=1:ne
        % Rotation
        % Eq. 10 of Fujita et al. (2006)
        R(:,:,i)=[cos(phase0(i)),-sin(phase0(i));sin(phase0(i)),cos(phase0(i))];
        
        % Scattering
        % Eq. 8 of Fujita et al. (2006)
        S(:,:,i)=[Sx(i) 0; 0 Sy(i)];
        % Rotate Scattering
        % Eq. 11 of Fujita et al. (2006)
        S(:,:,i)=R(:,:,i)*S(:,:,i)*R(:,:,i)';
        
        % Transmission
        % Eq. 5 of Fujita et al. (2006)
        T(:,:,i)=[T_ix(i) 0; 0 T_iy(i)];
        % Rotate transmission
        % Eq. 9 of Fujita et al. (2006)
        T(:,:,i)=R(:,:,i)*T(:,:,i)*R(:,:,i)';
    end
    
    % Calculate cumulative Transmission (product summation of Eq. 9)
    % cT(n)=Product{R' T R}i=1,n
    cT=T;
    for i=2:ne
        cT(:,:,i)=cT(:,:,i-1)*cT(:,:,i);
    end
    
    Rt=[cos(x(j)),sin(x(j));-sin(x(j)),cos(x(j))];
    Et=Rt*[1;0]; %transmiter in the frame of reference; [E_PT E_OT] 
    
    % Eq. 12 of Fujita et al. (2006)
    for i=1:ne
        Er(:,i)=cT(:,:,i)*S(:,:,i)*cT(:,:,i)*Et; %Receiver in the frame of reference
        Erant(:,i)=transpose(Rt)*Er(:,i); %Receiver in the frame of antennas
        Erp(i)=Erant(1,i);
        Ero(i)=Erant(2,i);
        
        % Or reconstruct signal for any azimuthal angle
        Err(:,:,i) = cT(:,:,i)*S(:,:,i)*cT(:,:,i)*Rt;
        Errant(:,:,i) = transpose(Rt)*Err(:,:,i);
        Shhn(i) = Errant(1,1,i);
        Svhn(i) = Errant(1,2,i);
        Shvn(i) = Errant(2,1,i);
        Svvn(i) = Errant(2,2,i); 
    end
    
    ErPar(:,j) = Erp;
    ErPer(:,j) = Ero;
    
    Shh(:,j) = Shhn;
    Svh(:,j) = Svhn;
    Shv(:,j) = Shvn;
    Svv(:,j) = Svvn;

end

%% Obtain power anomaly and phase difference

% Detrend power to get anomaly
PrPar=20*log10(abs(Shh));
PrPar=detrend(PrPar','constant')';
PrPar=AverageDepth(PrPar,z(2:end),dzPowerNorm);

PrPer=20*log10(abs(Shv));
PrPer=detrend(PrPer','constant')';
PrPer=AverageDepth(PrPer,z(2:end),dzPowerNorm);

% Get phase difference
for j=1:X-1
    diffphasepar(:,j)=angle(ErPar(:,j+1).*conj(ErPar(:,j)));
    diffphaseper(:,j)=angle(ErPer(:,j+1).*conj(ErPer(:,j)));
    thetaMid(j)=(x(j)+x(j+1))/2;
end

%% Calculate and plot polarimetric difference

% Calculate polarimetric coherence (hhvv)
C_hhvv = (Shh.*conj(Svv)) ./ (sqrt(abs(Shh).^2).*sqrt(abs(Svv).^2)); % hhvv coherence, Jordan et al. (2019) Eq. 15
C_hhvv = conj(C_hhvv); % De-ramp phase, Jordan et al. (2020) Eq. 7
P_hhvv = angle(C_hhvv); % hhvv phase difference, Jordan et al. (2019) Eq. 14
R = real(C_hhvv); I = imag(C_hhvv); 
[~,dRdz] = gradient(R,dx(1),dz(1)); 
[~,dIdz] = gradient(I,dx(1),dz(1)); 
Dz_hhvv = (R.*dIdz - I.*dRdz) ./ (R.^2 + I.^2); % Jordan et al. (2019) Eq. 23
P_hhvv_u = unwrap(P_hhvv,[],2); 

% Extract values from symmetry axis
for ii = 1:Z-1
    [~,phase0idx(ii)] = closest(x,phase0(ii));
end

% Back-calculate E2E1 (they should match perfectly)
for ii = 1:length(phase0idx)
    Dy_hhvv_phase0(ii) = Dz_hhvv(ii,phase0idx(ii));
end
CalculateE2E1 = @(P,ei,dei,c,fc) P .* (2*sqrt(ei)/dei) .* c/(4*pi*fc); 
E2E1 = CalculateE2E1(Dy_hhvv_phase0,ei,dei,c,fc);

%% Plot figure results

% Figure parameters
position = [0,0,1/2,1];
m = 1; % # subplot rows
n = 3; % # subplot columns
cLim1 = [-10 10];
cLim2 = [-6 +6];
cLim3 = [-pi +pi];
E1Color1 = [20, 90, 50]./255;
E1Color2 = [0, 255, 0] ./ 255;
E2Color1 = [125, 102, 8] ./ 255;
E2Color2 = [255, 255, 0] ./ 255;
xLim = [-90 90];
yLim = [0 1500];

% Set up figure structure 
% Colour scheme from https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
figA = figure; 
set(figA,'Units','Normalized','Position',position);
hold on
sp(1) = subplot(m,n,1); hold all, box on, grid on, cmocean('balance')
sp(2) = subplot(m,n,2); hold all, box on, grid on, cmocean('balance')
sp(3) = subplot(m,n,3); hold all, box on, grid on, cmocean('balance')
ax(1) = imagesc(sp(1),xdeg-91,z(2:end),PrPar);
ax(2) = imagesc(sp(2),thetaMid*180/pi-91,z(2:end),diffphasepar*180/pi);
ax(3) = imagesc(sp(3),xdeg-91,z(2:end),P_hhvv);

for ii = 1:3
    plot(sp(ii),phase0deg,z(2:end),'color',E1Color1,'lineWidth',2.5);
    plot(sp(ii),phase0deg-90,z(2:end),'color',E2Color1,'lineWidth',2.5);
    cb(ii) = colorbar(sp(ii),'SouthOutside');
    xlabel(sp(ii),['Azimuth (\theta) [' char(176) ']']);
end

% Set limits and ticks
xlim(sp,xLim)
ylim(sp,yLim)
set(sp,'Xtick',xLim(1):45:xLim(2))
set(sp,'XMinorTick','on')
for ii = 1:length(sp)
    sp(ii).XAxis.MinorTickValues = xLim(1):10:xLim(2);
end
cb(3).Ticks = cLim3(1):pi:cLim3(2);
cb(3).TickLabels = {'-\pi','0','+\pi'};
caxis(sp(1),cLim1);
caxis(sp(2),cLim2);
caxis(sp(3),cLim3);

set(sp,'YDir','reverse');
set(sp,'layer','top');

% Label axes
ylabel(sp(1),'Depth [m]')
ylabel(cb(1),'Relative power [dBm]');
ylabel(cb(2),['Angle [' char(176) ']']);
ylabel(cb(3),['{\it\phi_{hhvv}} [rad]']);
set(sp([2 3]),'YTickLabel',[])

% Assign titles
title(sp(1),{'{\bf(a)} Co-polarised power anomaly'},'FontWeight','Normal')
title(sp(2),{'{\bf(b)} Co-polarised phase difference'},'FontWeight','Normal')
title(sp(3),{'{\bf(c)} {\ithhvv} phase angle'},'FontWeight','Normal')

% Move subplot positions
tighten_subplots(sp(1:3),[m n],[0 0.03]) % 
for ii = 1:(m*n)  % Get original positions
    spPos(ii,:) = get(sp(ii),'Position'); 
    spPos(ii,3:4) = [0.2134 0.2584];
    set(sp(ii),'Position',spPos(ii,:))
end 

% Save output as figure
if saveData
    startup
    fileType = '.pdf';%'.pdf';
    MSdir = '2020Young_WAISD';
    folderOut = [gwd,'Manuscripts/',MSdir,'/'];
    fileOut = ['WAISD_fig3a',fileType];
    set(figA,'Color','w')
    export_fig(figA,[folderOut,fileOut],'-q101')
end

%% Function for averaging depth between two layers
function y=AverageDepth(x,z,dz)

for i=1:length(z)
    iok=find(z>z(i)-dz/2&z<z(i)+dz/2);
    y(i,:)=mean(x(iok,:),1);
end

end

%% Function to tighten subplots
function tighten_subplots(ax,mn,shift)
    
vShift = shift(1);
hShift = shift(2); 
hShift = -hShift; % To make vertical and horizontal directions consistent

% Transpose axes matrix if row-column order is reversed
if isa(ax,'matlab.graphics.axis.Axes') == 1
    if size(ax) ~= mn
        ax = ax';
        [m,n] = size(ax);
    end
end

m = mn(1); n = mn(2); 
oldspPosArray = get(ax,'Position');
oldspPosArray = reshape(oldspPosArray,m,n);
newspPosArray = cell(size(oldspPosArray,1),size(oldspPosArray,2)); 

% Move axes
for ii = 1:m
    for jj = 1:n
        oldspPos = oldspPosArray{ii,jj}; % Get original position of axis
        xn = oldspPos(1) + hShift*(jj-1); 
        yn = oldspPos(2)+vShift*(ii-1); 
        dxn = oldspPos(3); 
        dyn = oldspPos(4); 
        newspPos = [xn yn dxn dyn]; % Calculate new position of axis
        set(ax(ii,jj),'Position',newspPos) % Assign new position of axis
        newspPosArray(ii,jj) = {newspPos}; 
    end
end

end