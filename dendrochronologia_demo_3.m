% This is a demo script for the MATLAB version of the Vaganov-Shashkin model (VSM)
% This script produces the analyses for Yakutia/Indigirka in Anchukaitis et al. 2019 (submitted)
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

% get the Chokurdakh and remove the columns we won't use
chok   = textread('chokurdakh_daily_complete.txt'); chok(:,[1 5 6 7 9 10 11 13 14]) = [];
year   = chok(:,1); month = chok(:,2); day = chok(:,3);

doy = date2doy(year,month,day);

startYear = min(unique(year));
endYear   = max(unique(year));
allYears = [startYear:endYear]; 
nyears = length(allYears);

temperatureData   = NaN(366,nyears);
precipitationData = NaN(366,nyears);

for i = 1:nyears
    indx = find(year==allYears(i));
    temperatureData(doy(indx),i) = chok(indx,4);
    precipitationData(doy(indx),i) = chok(indx,5);
end

[T,P] = vsm_fillmiss(temperatureData,precipitationData,startYear,endYear);

T0 = T;
Tfil = vsm_filter(T0,'V'); % close mimic of the internal lowpass FORTRAN filter
phi = 70.0;

% call the parameter file
generic_parameters

%% Run the MATLAB mimic model
for j = 0:10  
    parameters.Tf(1) = j;
    output(j+1) = vsm(Tfil,P,phi,startYear,endYear,parameters);
    simulatedSeries(:,j+1) = output(j+1).trw';
end

simulatedSeries(simulatedSeries==0) = NaN;  % for plotting

% get the actual chronology file
load indigirka_chronology.mat

[commonYear,ia,ib] = intersect(allYears,Indigirka(:,1));
[R,P]   = corrcoef([simulatedSeries(ia,:) Indigirka(ib,2)],'rows','pairwise');
good    = find(commonYear~=1991);
[R1,P1] = corrcoef([simulatedSeries(ia(good),:) Indigirka(ib(good),2)],'rows','pairwise');


%% 
figure(1); clf
subplot(2,1,1)
ex1 = plot([startYear:endYear],standardize(simulatedSeries),'color',[1 0.8 0.8],'linewidth',1.5); hold on
ex2 = plot([startYear:endYear],standardize(simulatedSeries(:,2)),'r','linewidth',1.5); hold on
rx = plot(Indigirka(:,1),standardize(Indigirka(:,2)),'k','linewidth',1.5)
xlim([startYear-1 endYear+1])
xlabel('YEAR')
ylabel('INDEX')
ylim([-3.5 3])
set(gca,'XMinorTick','on','YMinorTick','on')
text(1945,2.4,'A','fontsize',14)
legend([ex1(1) ex2(1) rx],'ENSEMBLE','BEST','OBSERVED','location','southeast')
legend boxoff
title('LARCH, INDIGIRKA (YAKUTIA), RUSSIA')

subplot(2,2,3)
yyaxis left
ax1 = gca;
bx1 = bar([0:10],R(1:end-1,end));
bx1.BarWidth = 1;
bx1.FaceColor=[0 0 1];
ax1.YColor = [0.15 0.15 0.15];

ylabel('CORRELATION')
xlabel('T_{minimum}')

yyaxis right
ax2 = gca;
px1 = plot([0:10],P(1:end-1,end),'ro-','MarkerFaceColor','r'); hold on;
text(8,0.27,'P','color','r')
plot([0 10],[0.05 0.05],'r--')
set(gca,'YScale','log')
text(10.4,0.42,'B','fontsize',14)
ylabel('SIGNIFICANCE','rot',-90,'VerticalAlignment','bottom')
ax1.YColor = [1 0 0];

subplot(2,2,4)
bx2 = bar([0:10],100*((sum(isnan(simulatedSeries))-1)./length(simulatedSeries)));
bx2.BarWidth = 1;
bx2.FaceColor = [1 0.8 0.8];
set(gca,'YAxisLocation','right')
ylabel('% ZERO GROWTH YEARS','rot',-90,'VerticalAlignment','bottom')
xlabel('T_{minimum}')
text(-1,23,'C','fontsize',14)

print -depsc indigirka_panel.eps
