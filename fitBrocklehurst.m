% opts.format='pdf'; opts.outputDir='.'; publish('fitFreeman.m',opts);
clear all
close all

% load fig 2 data from 
%   Freeman, Michael L., et al. "The effect of pH on cell lethality induced by hyperthermic treatment." Cancer 45.9 (1980): 2291-2300.
fig2data = readtable('brocklehurstthermal.csv')

survival = min(1,max(1e-6,table2array(fig2data(strcmp(fig2data.runid,'survival') ,2:end))))
damage   = log(survival.^-1) + 1e-6
timedata    = table2array(fig2data(strcmp(fig2data.runid,'time') ,2:end)) * 60; % sec
temperature = table2array(fig2data(strcmp(fig2data.runid,'temperature') ,2:end)) +273 ; % K


% verify digitization
figure(1)
GasConst  = 8.314 ; % J/K / mol
plot((GasConst * temperature   ).^(-1),log(timedata.*damage.^-1),'x')

handlefit=figure(2)
% get Rsquared
mdl = fitlm((GasConst * temperature   ).^(-1),log(timedata.*damage.^-1) )
plot(mdl)
xlabel( '$\frac{1}{RT}$','interpreter','latex')
ylabel( '$\ln \frac{\Delta \tau}{\Omega}$','interpreter','latex')
title(sprintf('Arrhenius fit, R^2=%f, Ea=%9.2e, A=%9.2e',mdl.Rsquared.Ordinary, mdl.Coefficients.Estimate(2) , exp(abs( mdl.Coefficients.Estimate(1)))))
saveas(handlefit,'ArrheniusFitThermal','png')

