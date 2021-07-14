function fc = WAIS_firn_correction(xIn,eiPer)

% Function to correct for firn density, following Eq. A7 in Jordan et al. 
% (2020). 
%
% Note: Users will need to edit the inFolder and inFile variable to point
%   towards their copy of the WAIS Divide Firn Density dataset. This is
%   obtained from the NSIDC dataset of Gregory et al. (2014): 
%   https://nsidc.org/data/nsidc-0602/

% References
% Gregory et al. (2014) Impact of physical properties and accumulation rate 
%   on pore close-off in layered firn. The Cryosphere 8:91-105. 
%   doi:10.5194/tc-8-91-2014
% Jordan, T. M. et al (2018) Estimation of ice fabric within Whillans Ice
%   Stream using polarimetric phase-sensitive radar sounding. Annals of
%   Glaciology 61(81):74-83. doi:10.1017/aog.2020.6.
%
% TJ Young
% 01 July 2020

% Constants
Rice = 0.91675; % Density of ice [g cm^-3]
e = 2.71828; % Value of e

% Load data
%inFolder = '/Volumes/GoogleDrive/My Drive/ITGC-TIME/Datasets/WAIS Divide Firn Density/';
inFolder = 'H:/My Drive/ITGC-TIME/Datasets/WAIS Divide Firn Density/';
inFile = 'WAISDivide_FirnProperties.xlsx';
opts = detectImportOptions([inFolder,inFile]);
WF = readtable([inFolder,inFile],opts,'Sheet','WAIS Divide');
WF = WF(:,2:8);

if size(xIn,1) > size(xIn,2)
    xIn = xIn';
end

% Eyeball Fig. 9a data (somehow the data only goes between 55 - 75 m...)
addx = [0 10 20 30 40 50]';
addy = [0.44 0.52 0.59 0.65 0.7 0.74]';
addy = interp1(addx,addy,0:50)';
addx = [0:50]'; 

x = [addx;WF.DepthOfCenterOfSample];
y = [addy;WF.SampleDensity];

lnx = log(x);
lny = log(y);

[xData, yData] = prepareCurveData( lnx, lny );

% Set up fittype and options.
ft = fittype('poly2');
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';

% Fit model to data.
[fitresult,~] = fit( xData, yData, ft, opts );
xIn = log(xIn);
Rfirn = fitresult(xIn);

% Convert out of log
xIn = e.^xIn;
Rfirn = e.^Rfirn;
Rice = repmat(Rice,size(xIn))';

% Correct for log effects
if xIn(1) < e
    [~,cI] = closest(xIn,e^(1));
    if cI-1 > 4
        vq = interp1(xIn(cI:end),Rfirn(cI:end),xIn(1:cI),'pchip','extrap');
    elseif (cI - 1 >= 1) && (cI - 1 < 5)
        vq = interp1(xIn(cI:end),Rfirn(cI:end),xIn(1:cI),'linear','extrap');
    elseif cI == 1
        vq = Rfirn(2) - diff(Rfirn(2:3));
    end
    Rfirn(1:cI,:) = vq';
end

% Convert to firn correction parameter (Eq. A7 in Jordan et al. 2020)
v = Rfirn ./ Rice; % Convert density to ice volume fraction
v(v>1) = 1; % Set upper limit of v to 1
fv = @(v,Eper) (1/Eper) .* ( (v.^3.*Eper) + (2/3.*v.^2.*(1-v).*(Eper.^(2/3))) + (1/3.*v.*(1-v).^2.*(Eper.^(1/3))) );
fc = fv(v,eiPer);

if size(fc,1) > size(fc,2)
    v = v';
    fc = fc';
end

