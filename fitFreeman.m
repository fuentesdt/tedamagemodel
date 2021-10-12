clear all
close all


fig2data = readtable('freemanFig2.csv')

base10min = fig2data(strcmp(fig2data.group,'10minpH7.5') ,:);
acid10min = fig2data(strcmp(fig2data.group,'10minpH6.65'),:);
base15min = fig2data(strcmp(fig2data.group,'15minpH7.5') ,:);
acid15min = fig2data(strcmp(fig2data.group,'15minpH6.65'),:);
base20min = fig2data(strcmp(fig2data.group,'20minpH7.5') ,:);
acid20min = fig2data(strcmp(fig2data.group,'20minpH6.65'),:);


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


disp('build objective function')
mycostfcn = log( base10min.survival.^(-1)) - ...
   deltaTheat10min * frequencyfactor * exp(-(Ea0+Ea1* pHheat)/(GasConst * Theat10min)) - ...
   base10min.minute * frequencyfactor * exp(-(Ea0+Ea1* pHincbase)/(GasConst * Tinc))


disp('create optim prob')
convprob = optimproblem('Objective',mycostfcn );
%show(convprob)
problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
%% 
% Solve the new problem. The solution is essentially the same as before.
myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed');
[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'solver','lsqnonlin' )
