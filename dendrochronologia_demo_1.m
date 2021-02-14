% This is a demo script for the MATLAB version of the Vaganov-Shashkin model (VSM)
% This script produces the analyses for the White Mountains in Anchukaitis et al. 2019 (submitted)
clear; close all; clc;  % clear the workspace

myColormap = [... % 'Reds' from ColorBrewer
    0.9882    0.7333    0.6314
    0.9882    0.5725    0.4471
    0.9843    0.4157    0.2902
    0.9373    0.2314    0.1725
    0.7961    0.0941    0.1137
    0.6471    0.0588    0.0824
    0.4039         0    0.0510];

% set some plotting defaults
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','normal')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultTextFontWeight','normal') 
set(0,'defaultaxesfontsize',12); set(0,'defaulttextfontsize',12); 

% load the actual WM1 Chronology
load('wm1_actual.mat'); % Methuselah Walk

% call the script that sets up the parameter structure
wm1_parameter_setup

% get the input data
load('049632_input.mat');

% adjust temperatures for difference in site vs. meteorological station data
T  = T + 2.6; % estimated for station vs site elevation difference and adiabatic lapse rate
T0 = T;

% Mimic the 'V' filter from fortran (LOW-PASS HALF FILTER)
Tfil = vsm_filter(T0,'V'); % mimic of the lowpass FORTRAN filter used in the FORTRAN version

% Which temperature input series to use?
T = Tfil; % use the lowpass filtered temperature

%% Run the MATLAB mimic model
tic % based on testing, this loop may take anywhere from 1 to 20 seconds, depending on your system
output = vsm(T,P,phi,syear,eyear,parameters);
toc %  

% calculate correlation between actual chronology and simulated
[Ro] = corrcoef([chronology(:,2) output.trw']);

%% diagnostic figures
figure(1); clf
subplot(2,1,1)
plot([syear:eyear],zscore(output.trw),'r','linewidth',1.5); hold on
plot(chronology(:,1),zscore(chronology(:,2)),'k','linewidth',1.5)
legend('SIMULATED','OBSERVED','location','southeast')
legend boxoff
title('WHITE MOUNTAINS, METHUSELAH WALK')
xlabel('YEAR')
ylabel('INDEX')
set(gca,'xminortick','on','yminortick','on')
text(1955.5,1.5,'A','fontsize',14)
text(1957,1.2,['r = ',num2str(Ro(1,2),'%1.2f') ', p < 0.001'],'fontsize',12) %  num2str(Po(1,2),'%1.5f

subplot(2,2,4)
plot(mean(output.Gr,2),'k'); hold on;
plot(mean(output.GrW,2),'b'); 
plot(mean(output.GrT,2),'r'); 
lx2 = legend('G','g_W','g_T','location','west'); legend boxoff; % 
xlim([1 365]); ylim([0 1])
title('GROWTH RATES')
xlabel('DAY OF YEAR')
ylabel('GROWTH RATE','rot',-90,'VerticalAlignment','bottom')
set(gca,'yaxislocation','right')
set(gca,'xminortick','on','yminortick','on')
text(15,0.9,'C','fontsize',14)

subplot(2,2,3)
snx = plot(mean(output.snowdepth/100,2),'k'); hold on;
smx = plot(mean(output.sm,2),'b'); 
tnx = plot(mean(output.transpiration,2),'r'); 
[lx,icons]  = legend([snx tnx smx],'SNOW/100 (mm)','TRANS (mm/d)','SM (v/v)','location','northwest','AutoUpdate','off'); legend boxoff;
set(lx,'FontSize',8)
pos = get(lx,'Position'); set(lx,'Position',[pos(1)-0.015 pos(2)+0.015 pos(3) pos(4)])
xlim([1 365]); 
title('ENVIRONMENT')
xlabel('DAY OF YEAR')
set(gca,'xminortick','on','yminortick','on')
text(335,2.75,'B','fontsize',14)

print -depsc white_panel1.eps

clear Ro Po

%% test sensitivity to soil drainage rate
% you may wish to use a coarser valueStep at first in order to reduce the time spent in the loop
valueStep = 0.0001;
drainageRates = 0.001:valueStep:0.0200; 

tic  % based on testing, this may take anywhere from 20 to 120 seconds
outputm = NaN(length(syear:eyear),size(drainageRates,2));
for i = 1:length(drainageRates)
    parameters.rated = drainageRates(i);
    output(i) = vsm(T,P,phi,syear,eyear,parameters);
    [Ro(1:2,1:2,i),Po(1:2,1:2,i)] = corrcoef([chronology(:,2) output(i).trw']);
    outputm(:,i) = output(i).trw';
end
toc 

R0 = squeeze(Ro(1,2,:));
bestSimulation = find(R0==max(R0)); bestSimulation = bestSimulation(1);

%% 
figure(2); clf
subplot(2,1,1)
ex1 = plot([syear:eyear],zscore(outputm),'color',[1 0.8 0.8],'linewidth',1); hold on
ex2 = plot([syear:eyear],zscore(outputm(:,bestSimulation)),'color',[1 0 0],'linewidth',1.5); hold on
rx = plot(chronology(:,1),zscore(chronology(:,2)),'k','linewidth',1.5);
legend([ex1(1) ex2(1) rx],'ENSEMBLE','BEST','OBSERVED','location','southeast')
legend boxoff
title('WHITE MOUNTAINS, METHUSELAH WALK')
xlabel('YEAR')
ylabel('INDEX')
set(gca,'xminortick','on','yminortick','on')
text(1955.5,2.5,'A','fontsize',14)

subplot(2,2,3)
scatter(drainageRates,squeeze(Ro(1,2,:)),10,squeeze(Po(1,2,:)),'fill')
caxis([0 0.05])
cbx = colorbar('location','East','YDir', 'reverse');
colormap(flipud(myColormap))
xlim([0 0.021])
xlabel('DRAINAGE RATE')
ylabel('CORRELATION')
pos = get(cbx,'position');
set(cbx,'position',[pos(1)+0.06 pos(2) 0.03 pos(4)])
title(cbx,'p')
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 0.28 pos(4)])
text(0.01875,0.75,'B','fontsize',14)
set(gca,'xminortick','on','yminortick','on','box','on')

subplot(2,2,4)
hist(squeeze(Ro(1,2,:)))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w')
xlabel('CORRELATION')
ylabel('FREQUENCY','rot',-90,'VerticalAlignment','bottom')
set(gca,'yaxislocation','right')
xlim([0.0 1])
text(.05,45,'C','fontsize',14)
set(gca,'xminortick','on','yminortick','on')

print -depsc white_panel2.eps

