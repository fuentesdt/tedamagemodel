% opts.format='pdf'; opts.outputDir='.'; publish('fitFreeman.m',opts);
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
logA = optimvar('logA','LowerBound',0);
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
   deltaTheat10min  * exp(logA -(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   base10min.minute * exp(logA -(Ea0+Ea1* pHincbase)/(GasConst * Tinc))).^2) + ...
                 sum((log( acid10min.survival.^(-1)) - ...
   deltaTheat10min  * exp(logA -(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   acid10min.minute * exp(logA -(Ea0+Ea1* pHincacid)/(GasConst * Tinc))).^2) ;
mycostfcn15min = sum((log( base15min.survival.^(-1)) - ...
   deltaTheat15min  * exp(logA -(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   base15min.minute * exp(logA -(Ea0+Ea1* pHincbase)/(GasConst * Tinc))).^2) + ...
                 sum((log( acid15min.survival.^(-1)) - ...
   deltaTheat15min  * exp(logA -(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   acid15min.minute * exp(logA -(Ea0+Ea1* pHincacid)/(GasConst * Tinc))).^2) ;
mycostfcn20min = sum((log( base20min.survival.^(-1)) - ...
   deltaTheat20min  * exp(logA -(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   base20min.minute * exp(logA -(Ea0+Ea1* pHincbase)/(GasConst * Tinc))).^2) + ...
                 sum((log( acid20min.survival.^(-1)) - ...
   deltaTheat20min  * exp(logA -(Ea0+Ea1* pHheat)/(GasConst * Theat     )) - ...
   acid20min.minute * exp(logA -(Ea0+Ea1* pHincacid)/(GasConst * Tinc))).^2) ;
mycostfcn = mycostfcn10min + mycostfcn15min + mycostfcn20min;
show(mycostfcn )


disp('create optim prob')
convprob = optimproblem('Objective',mycostfcn );
%show(convprob)
problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
%% 
% Solve the new problem. The solution is essentially the same as before.
%'Algorithm', 'trust-region-reflective';
myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed','OptimalityTolerance',1.e-12,'FunctionTolerance', 1.000000e-12);
%myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed','OptimalityTolerance',1.e-12,'FunctionTolerance', 1.000000e-12,'Algorithm' , 'levenberg-marquardt');


x0.Ea0 = 6.28e5; % J/mol
x0.Ea1 = 0;
x0.logA = log(3.1e98);
[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'solver','lsqnonlin' )

initobjfunc = evaluate(mycostfcn ,x0)
optobjefunc = evaluate(mycostfcn ,popt)

%evaluate fit
damagepredictionbase10min = deltaTheat10min * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - base10min.minute * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHincbase)/(GasConst * Tinc));
damagepredictionacid10min = deltaTheat10min * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - acid10min.minute * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHincacid)/(GasConst * Tinc));
damagepredictionbase15min = deltaTheat15min * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - base15min.minute * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHincbase)/(GasConst * Tinc));
damagepredictionacid15min = deltaTheat15min * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - acid15min.minute * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHincacid)/(GasConst * Tinc));
damagepredictionbase20min = deltaTheat20min * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - base20min.minute * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHincbase)/(GasConst * Tinc));
damagepredictionacid20min = deltaTheat20min * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHheat)/(GasConst * Theat     )) - acid20min.minute * exp(popt.logA ) * exp(-(popt.Ea0+popt.Ea1* pHincacid)/(GasConst * Tinc));


figure(2)
plot(log( base10min.survival.^(-1)), damagepredictionbase10min , 'x',log( acid10min.survival.^(-1)), damagepredictionacid10min , 'x', log( base15min.survival.^(-1)), damagepredictionbase15min , 'x',log( acid15min.survival.^(-1)), damagepredictionacid15min , 'x', log( base20min.survival.^(-1)), damagepredictionbase20min , 'x',log( acid20min.survival.^(-1)), damagepredictionacid20min , 'x' )
xlabel( 'measured damage')
ylabel( 'predicted damage')

measureddamage =  [log( base10min.survival.^(-1));
                   log( acid10min.survival.^(-1));
                   log( base15min.survival.^(-1));
                   log( acid15min.survival.^(-1));
                   log( base20min.survival.^(-1));
                   log( acid20min.survival.^(-1))];

predicteddamage =[damagepredictionbase10min;
                  damagepredictionacid10min;
                  damagepredictionbase15min;
                  damagepredictionacid15min;
                  damagepredictionbase20min;
                  damagepredictionacid20min];

% get Rsquared
mdl = fitlm(measureddamage ,predicteddamage )
handlefit=figure(3)
plot(mdl)
xlabel( 'measured damage')
ylabel( 'predicted damage')
title(sprintf('Arrhenius fit, R^2=%f, Ea0=%9.2e, Ea1=%9.2e',mdl.Rsquared.Ordinary,popt.Ea0,popt.Ea1))
saveas(handlefit,'ArrheniusFit','png')
