% opts.format='pdf'; opts.outputDir='.'; publish('tcadamage.m',opts);
clear all
close all


% load fig 2 data from 
%   Freeman, Michael L., et al. "The effect of pH on cell lethality induced by hyperthermic treatment." Cancer 45.9 (1980): 2291-2300.
fig2data = readtable('freemanFig2.csv')

base10min = fig2data(strcmp(fig2data.group,'10minpH7.5') ,:);
acid10min = fig2data(strcmp(fig2data.group,'10minpH6.65'),:);
base15min = fig2data(strcmp(fig2data.group,'15minpH7.5') ,:);
acid15min = fig2data(strcmp(fig2data.group,'15minpH6.65'),:);
base20min = fig2data(strcmp(fig2data.group,'20minpH7.5') ,:);
acid20min = fig2data(strcmp(fig2data.group,'20minpH6.65'),:);


% verify digitization
figure(1)
semilogy(base10min.minute,base10min.survival)
hold
semilogy(acid10min.minute,acid10min.survival)
semilogy(base15min.minute,base15min.survival)
semilogy(acid15min.minute,acid15min.survival)
semilogy(base20min.minute,base20min.survival)
semilogy(acid20min.minute,acid20min.survival)


% setup curve fit
Ea0 = optimvar('Ea0','LowerBound',0);
Ea1 = optimvar('Ea1','LowerBound',0);
frequencyfactor = 3.1e98;
deltaTheat10min = [10;10;10];
deltaTheat15min = [15;15;15];
deltaTheat20min = [20;20;20];
pHheat = 7.5;
pHincbase = 7.5;
pHincacid = 6.65;
GasConst  = 8.314 ; % J/K / mol
Theat      = 45.5 + 273; % K
Tinc       = 37.0 + 273; % K

disp('build objective function')
mycostfcn10min = sum((log( base10min.survival.^(-1)) - ...
   deltaTheat10min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   base10min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincbase)/(GasConst * Tinc))).^2) + ...
                 sum((log( acid10min.survival.^(-1)) - ...
   deltaTheat10min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   acid10min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincacid)/(GasConst * Tinc))).^2) ;
mycostfcn15min = sum((log( base15min.survival.^(-1)) - ...
   deltaTheat15min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   base15min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincbase)/(GasConst * Tinc))).^2) + ...
                 sum((log( acid15min.survival.^(-1)) - ...
   deltaTheat15min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   acid15min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincacid)/(GasConst * Tinc))).^2) ;
mycostfcn20min = sum((log( base20min.survival.^(-1)) - ...
   deltaTheat20min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   base20min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincbase)/(GasConst * Tinc))).^2) + ...
                 sum((log( acid20min.survival.^(-1)) - ...
   deltaTheat20min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   acid20min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincacid)/(GasConst * Tinc))).^2) ;
mycostfcn = mycostfcn10min + mycostfcn15min + mycostfcn20min;
show(mycostfcn )


disp('create optim prob')
convprob = optimproblem('Objective',mycostfcn );
%show(convprob)
problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
%% 
% Solve the new problem. The solution is essentially the same as before.
myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed');
x0.Ea0 = 6.28e5; % J/mol
x0.Ea1 = 0;
[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'solver','lsqnonlin' )


%evaluate fit
survivalpredictionbase10min = deltaTheat10min * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - base10min.minute * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHincbase)/(GasConst * Tinc))
survivalpredictionacid10min = deltaTheat10min * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - acid10min.minute * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHincacid)/(GasConst * Tinc))
survivalpredictionbase15min = deltaTheat15min * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - base15min.minute * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHincbase)/(GasConst * Tinc))
survivalpredictionacid15min = deltaTheat15min * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - acid15min.minute * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHincacid)/(GasConst * Tinc))
survivalpredictionbase20min = deltaTheat20min * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - base20min.minute * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHincbase)/(GasConst * Tinc))
survivalpredictionacid20min = deltaTheat20min * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - acid20min.minute * frequencyfactor * exp(-(popt.Ea0+popt.Ea1* pHincacid)/(GasConst * Tinc))


figure(2)
plot(log( base10min.survival.^(-1)), survivalpredictionbase10min , 'x',log( acid10min.survival.^(-1)), survivalpredictionacid10min , 'x', log( base15min.survival.^(-1)), survivalpredictionbase15min , 'x',log( acid15min.survival.^(-1)), survivalpredictionacid15min , 'x', log( base20min.survival.^(-1)), survivalpredictionbase20min , 'x',log( acid20min.survival.^(-1)), survivalpredictionacid20min , 'x' )
