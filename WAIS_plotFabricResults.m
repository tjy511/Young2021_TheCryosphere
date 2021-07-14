function F = WAIS_plotFabricResults(E,S,HHVV,X,A,cfg)
% Plots results in a tiled format. This script requires the cmocean 
% function, which can be obtained from Thyng et al. (2016).
% 
% References
% Thyng, K. M. et al. (2016) True colors of oceanography: Guidelines for 
%   effective and accurate colormap selection. Oceanography 29(3):9-13. 
%   doi:10.5670/oceanog.2016.66
% 
% Input: 
% - E: Output of WAIS_ScatteringMatrix.
% - S: Output of WAIS_ReconstructAzimuth. 
% - HHVV: Output of WAIS_CalculateHHVV. 
% - X: Output of WAIS_FabricOrientation. 
% - A: Output of WAIS_FabricStrength. 
% - cfg: File containing all configuration parameters.
%
% Output: 
% - F: Figure handle. 
% 
% TJ Young
% 02 April 2021

%% Extract relevant parameters from structures
Azimuth = E.Azimuth;
Depth = E.Depth;
Azimuth_HHVV = HHVV.Azimuth; 
Depth_HHVV = HHVV.Depth;
PNPar = S.PNPar; 
PNPer = S.PNPer; 
PDPar = S.PDPar; 
C_HHVV = HHVV.C; 
P_HHVV = HHVV.P;
dP_HHVV = HHVV.dP;

% Define figure parameters
m = 2; 
n = 4; 
xLim = cfg.AzimuthRange;
yLim = cfg.rangeLim;

% Circ-shift data if domain is not [0 180]
cs = @(c,Dx) (circshift(c(1:end-1,:),Dx));
if cfg.AzimuthRange(1) ~= 0
    dx = mean(diff(Azimuth));
    Dx = (xLim(1)-Azimuth(1))/dx;
    Azimuth = xLim(1):dx:xLim(2);
    PNPar = cs(PNPar,Dx); PNPar = [PNPar;PNPar(1,:)];
    PNPer = cs(PNPer,Dx); PNPer = [PNPer;PNPer(1,:)];
    PDPar = cs(PDPar,Dx); PDPar = [PDPar;PDPar(1,:)];
    dx = mean(diff(Azimuth_HHVV));
    Dx = (xLim(1)-Azimuth_HHVV(1))/dx;
    Azimuth_HHVV = xLim(1):dx:xLim(2); 
    C_HHVV = cs(C_HHVV,Dx); C_HHVV = [C_HHVV;C_HHVV(1,:)];
    P_HHVV = cs(P_HHVV,Dx); P_HHVV = [P_HHVV;P_HHVV(1,:)];
    dP_HHVV = cs(dP_HHVV,Dx); dP_HHVV = [dP_HHVV;dP_HHVV(1,:)];
end

% Construct figure architecture and set figure params
F = figure; hold on
for ii = 1:m
    for jj = 1:n
        sp((ii-1)*n + jj) = subplot(m,n,(ii-1)*n + jj); 
        hold on, box on, grid on
    end
end
set(F,'position',[10+690,10,700,700]);
set(sp,'layer','top');

% (a) A-scope
ii = 1;
h1 = plot(sp(ii),db(abs(E.Ehh)),Depth);
h2 = plot(sp(ii),db(abs(E.Evh)),Depth);
h3 = plot(sp(ii),db(abs(E.Ehv)),Depth);
h4 = plot(sp(ii),db(abs(E.Evv)),Depth);
xlabel(sp(ii),'Power [dB (V_{RMS})]')
legend(sp(ii),[{'\its_{hh}','\its_{vh}','\its_{hv}','\its_{vv}'}],'location','southeast');

% (b) Co-parallel (HH) power anomaly
ii = 2; 
h = imagesc(sp(ii),Azimuth,Depth,PNPar');
cmocean(sp(ii),'balance')
caxis(sp(ii),[-5 +5])
xlim(sp(ii),xLim)

% (c) Cross-parallel (VH) power anomaly
ii = 3;
h = imagesc(sp(ii),Azimuth,Depth,PNPer');
cmocean(sp(ii),'balance')
caxis(sp(ii),[-5 +5])
xlim(sp(ii),xLim)

% (d) Co-parallel (HH) phase difference
ii = 4;
h = imagesc(sp(ii),Azimuth,Depth,rad2deg(PDPar)');
cmocean(sp(ii),'balance')
caxis(sp(ii),[-3 +3])
xlim(sp(ii),xLim)

% (e) HHVV coherence
ii = 5;
h = imagesc(sp(ii),Azimuth_HHVV,Depth_HHVV,abs(C_HHVV')); 
cmocean(sp(ii),'-grey')
caxis(sp(ii),[0 1])

% (f) HHVV Phase angle
ii = 6;
h = imagesc(sp(ii),Azimuth_HHVV,Depth_HHVV,P_HHVV'); 
cmocean(sp(ii),'balance')
caxis(sp(ii),[-pi +pi]);

% (g) HHVV Phase derivative
ii = 7;
h = imagesc(sp(ii),Azimuth_HHVV,Depth_HHVV,dP_HHVV'); 
cmocean(sp(ii),'balance')
caxis(sp(ii),[-0.1 +0.1]);

% (h) Azimuthal fabric strength
ii = 8;
h = scatter(sp(ii),A.a,Depth_HHVV,10,'k','filled'); 
%k = erbar(sp(ii),A.a,Depth_HHVV,+A.a_se,-A.a_se);
k = errorbar(sp(ii),A.a,Depth_HHVV,A.a_se,...
    'horizontal','.','MarkerSize',1,'Color','k','LineWidth',0.5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k',...
    'Color',[1 1 1].*0.7,'CapSize',0);
xlim(sp(ii),[0 0.5]);
xlabel(sp(ii),'\itE_2 - E_1')

% Plot eigenvectors
alpha = 0.25;
for ii = 2:7
    h = scatter(sp(ii),wrapTo180(X.E1.*2)./2,Depth,1,'g','filled',...
        'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha); % E1
    k = scatter(sp(ii),wrapTo180(X.E2.*2)./2,Depth,1,'y','filled',...
        'MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha); % E2
end

% Adjust axes
for ii = 1:8
    set(sp(ii),'YLim',yLim);
end

for ii = 2:7
    cb(ii) = colorbar(sp(ii),'southoutside');
    sp(ii).XAxis.MinorTickValues = xLim(1):10:xLim(2);
    set(sp(ii),'XMinorTick','on')
    set(sp(ii),'XTick',xLim(1):45:xLim(2))
    set(sp(ii),'XLim',xLim,'YLim',yLim);
    xlabel(sp(ii),['Azimuth (\theta) [' char(176) ']']);
end

% Set colour bar
cb(4).Ticks = -3:3:3;
cb(6).Ticks = -pi:pi:+pi;
cb(6).TickLabels = {'-\pi','0','+\pi'};
ylabel(cb(2),'Relative power [dBm]')
ylabel(cb(3),'Relative power [dBm]')
ylabel(cb(4),'Angle [rad]')
ylabel(cb(5),'|{\itc_{hhvv}}|')
ylabel(cb(6),'\phi_{\ithhvv} [rad]')
ylabel(cb(7),'{\itd}\phi_{\ithhvv}/{\itdz} [rad m^{-1}]')

ylabel(sp(1),'Depth [m]'); ylabel(sp(5),'Depth [m]');
set(sp([2 3 4 6 7 8]),'YTickLabel',[])
set(sp,'YDir','reverse');

% Label subplots
title(sp(1),{'{\bf(a)} A-scope';' '},'FontWeight','Normal')
title(sp(2),{'{\bf(b)} HH power anomaly';' '},'FontWeight','Normal')
title(sp(3),{'{\bf(c)} VH power anomaly';' '},'FontWeight','Normal')
title(sp(4),{'{\bf(d)} HH phase difference';' '},'FontWeight','Normal')
title(sp(5),{'{\bf(e)} {\ithhvv} coherence';' '},'FontWeight','Normal')
title(sp(6),{'{\bf(f)} {\ithhvv} phase angle';' '},'FontWeight','Normal')
title(sp(7),{'{\bf(g)} {\ithhvv} phase gradient';' '},'FontWeight','Normal')
title(sp(8),{'{\bf(h)} Azimuthal';'fabric strength';},'FontWeight','Normal')

% Adjust subplot size and positions
tighten_subplots([sp(1:4);sp(5:8)],[m n],[0.06 0.02])
spPos(1,:) = get(sp(1),'Position'); spPos(2,:) = get(sp(2),'Position'); 
set(sp(1),'Position',[spPos(1,1) spPos(2,2) spPos(1,3) spPos(2,4)])
spPos(7,:) = get(sp(7),'Position'); spPos(8,:) = get(sp(8),'Position'); 
set(sp(8),'Position',[spPos(8,1) spPos(7,2) spPos(8,3) spPos(7,4)])

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
