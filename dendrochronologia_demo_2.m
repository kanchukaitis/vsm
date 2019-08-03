% This is a demo script for the MATLAB version of the Vaganov-Shashkin model (VSM)
clear; close all; clc; % clean up the workspace

myColormap = [... % 'Reds' from ColorBrewer
    1.0000    0.9608    0.9412
    0.9961    0.8784    0.8235
    0.9882    0.7333    0.6314
    0.9882    0.5725    0.4471
    0.9843    0.4157    0.2902
    0.9373    0.2314    0.1725
    0.7961    0.0941    0.1137
    0.6471    0.0588    0.0824
    0.4039         0    0.0510];

%% read in the climate data
% the file will have to be in your current directory or your path
% in the workspace, make sure you have P, T, syear, eyear, and phi
load('305426_mohonk.mat');

% get the Mohonk hemlock chronology, NY004r, Cook and Jacoby 1977
mohonk

% smooth the input temperature to follow Vaganov et al. 2011
T = vsm_filter(T,'A');

%% read in the parameter file
a06_parameters % parameter file from Anchukaitis et al. 2006
parameters.rated = 0.0115;  % change drainage rate for the 'Humpty Dumpty' rocky slope at Mohonk - Vaganov et al. 2011

% use a Latin Hypercube design to sample from the parameter space for Tf(1)
% and Tf(2)
X = lhsdesignbnd(1000,2,[0 11],[10 20],[false false]);

% overlapping period of meteorological data and tree-ring chronology
[~,idx1,idx2] = intersect(crn(:,1),syear:eyear);

for i = 1:length(X)
     % [i]
     parameters.Tf(1) = X(i,1);
     parameters.Tf(2) = X(i,2);
     output(i) = vsm(T,P,phi,syear,eyear,parameters);
     [Ro(1:2,1:2,i),Po(1:2,1:2,i)] = corrcoef([crn(idx1,2) output(i).trw(idx2)']);
     outputt2(:,i) = output(i).trw';
end

% load mohonk_temperature_ensemble_lh2.mat

for j = 1:1000
    meanGrowthRate(:,j) = nanmean(output(j).Gr,2);
    meanGrowthRateT(:,j) = nanmean(output(j).GrT,2);
    meanGrowthRateW(:,j) = nanmean(output(j).GrW,2); 
end

R0 = squeeze(Ro(1,2,:))
bestSimulation = find(R0==max(R0)); bestSimulation = bestSimulation(1);

%%
figure(2); clf;
subplot(2,1,1)
ex2 = plot([syear:eyear],zscore(outputt2),'color',[1 0.7 0.7],'linewidth',1); hold on
ex = plot([syear:eyear],zscore(outputt2(:,bestSimulation)),'color',[1 0 0],'linewidth',1.5); hold on
rx = plot(crn(:,1),zscore(crn(:,2)),'k','linewidth',1.5);
lx = legend([ex2(1) ex(1) rx],'ENSEMBLE','BEST','REAL','location','southwest')
legend boxoff
set(lx,'Position',[0.22 0.81 0.18 0.10])
xlim([1890 2005])
title('EASTERN HEMLOCK, MOHONK, SHAWANGUNK MOUNTAINS')
xlabel('YEAR')
set(gca,'xminortick','on','yminortick','on')
text(1891,3.5,'A','fontsize',14)

subplot(2,2,3)
scatter(X(:,1),X(:,2),50,R0,'filled')
caxis([0.3 0.6])
cbx = colorbar('location','East');
colormap(myColormap)
xlim([0 10])
ylim([11 20])
xlabel('T_{minimum}')
ylabel('T_{lower optimum}')
pos = get(cbx,'position');
set(cbx,'position',[pos(1)+0.06 pos(2) 0.03 pos(4)])
title(cbx,'r')
pos = get(gca,'position');
set(gca,'position',[pos(1) pos(2) 0.28 pos(4)])
text(0.5,19,'B','fontsize',14) %,'BackgroundColor','w')
set(gca,'xminortick','on','yminortick','on','box','on')
title('PARAMETER SAMPLING')

subplot(2,2,4)
px1 = plot(meanGrowthRate,'color',[0.8 0.8 0.8]); hold on;
px2 = plot(meanGrowthRateT,'color',[1 0.8 0.8],'linewidth',0.5);
px3 = plot(meanGrowthRateT(:,bestSimulation),'color',[1 0 0],'linewidth',1.5);
px5 = plot(meanGrowthRateW(:,bestSimulation),'color',[0 0 1],'linewidth',1.5);
px4 = plot(meanGrowthRate(:,bestSimulation),'k','linewidth',1.5); hold on;
tickday =  [1; 31; 59; 90; 120; 151; 181; 212; 243; 273; 304; 334; 365];
set(gca,'XTick',tickday(1:12)+15)
set(gca,'XTickLabel',{'JAN'; ' '; 'MAR'; ' '; 'MAY'; ' '; 'JUL '; ' '; 'SEPT'; '  '; '  NOV'; ' '},'FontSize',10)
lx2 = legend([px4(1) px3(1) px5(1)],'G','g_T','g_W','location','south'); legend boxoff; 
xlim([1 365]); ylim([0 1])
title('GROWTH RATES','FontSize',12)
xlabel('DAY OF YEAR','FontSize',12)
set(gca,'xminortick','on','yminortick','on')
text(15,0.9,'C','fontsize',14)

% save mohonk_temperature_ensemble_lh2.mat
printeps('mohonk_panel.eps',2)
