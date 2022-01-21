% opts.format='pdf'; opts.outputDir='.'; publish('fitFreeman.m',opts);
clear all
close all

% load 2 data 
summarydata = readtable('summarydata.csv');
summarydata.Temperature = summarydata.Temperature + 273;
%summarydata.pH = summarydata.pH ;
summarydata.Time = summarydata.Time * 60 ;
summarydata.Viabilty = summarydata.Viabilty +1.e-6 ;

temperaturesubset = summarydata(strcmp(summarydata.expid,'TempOnly'),:)
pHsubset = summarydata(strcmp(summarydata.expid,'pHOnly'),:)



% verify digitization
figure(1)
plot(summarydata.Time,summarydata.pH,'b')
figure(2)
plot(summarydata.Time,summarydata.Temperature,'k')
figure(3)
plot(summarydata.Time,summarydata.Viabilty,'g')


% setup curve fit
Ea0 = optimvar('Ea0','LowerBound',0);
Ea1 = optimvar('Ea1','LowerBound',0);
logA = optimvar('logA','LowerBound',0);
GasConst  = 8.314 ; % J/K / mol

disp('build objective function')
fullcostfcn = sum((log( summarydata.Viabilty.^(-1)) - ...
   summarydata.Time.* exp(logA -Ea0*(GasConst * summarydata.Temperature +Ea1* summarydata.pH ).^(-1)) ).^2)  
phcostfcn = sum((log( pHsubset.Viabilty.^(-1)) - ...
   pHsubset.Time.* exp(logA -Ea0*(GasConst * pHsubset.Temperature +Ea1* pHsubset.pH ).^(-1)) ).^2)  
tempcostfcn = sum((log( temperaturesubset.Viabilty.^(-1)) - ...
   temperaturesubset.Time.* exp(logA -(Ea0)./(GasConst * temperaturesubset.Temperature  )) ).^2)  

disp('create optim prob')
convprob = optimproblem('Objective',fullcostfcn );
phprob   = optimproblem('Objective',  phcostfcn );
tempprob = optimproblem('Objective',tempcostfcn );
%show(convprob)
%% problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjective');
%% 
% Solve the new problem. The solution is essentially the same as before.
%'Algorithm', 'trust-region-reflective';
myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed','OptimalityTolerance',1.e-12,'FunctionTolerance', 1.000000e-12);
%myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed','OptimalityTolerance',1.e-12,'FunctionTolerance', 1.000000e-12,'Algorithm' , 'levenberg-marquardt');


x0.Ea0 = 6.28e5; % J/mol
x0.Ea1 = 0;
x0.logA = log(3.1e98);
[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'solver','lsqnonlin' )
[poptph,fvalph,exitflagph,outputph] = solve(phprob,x0,'Options',myoptions, 'solver','lsqnonlin' )
[popttemp,fvaltemp,exitflagtemp,outputtemp] = solve(tempprob,x0,'Options',myoptions, 'solver','lsqnonlin' )

initobjfunc = evaluate(fullcostfcn ,x0)
optobjefunc = evaluate(fullcostfcn ,popt)

%evaluate fit
damageprediction= exp(popt.logA ) * summarydata.Time.*exp(-popt.Ea0*(GasConst * summarydata.Temperature+popt.Ea1* summarydata.pH).^(-1));
measureddamage =  log( summarydata.Viabilty.^(-1));

damagepredictionph= exp(poptph.logA ) * pHsubset.Time.*exp(-poptph.Ea0*(GasConst * pHsubset.Temperature+poptph.Ea1* pHsubset.pH).^(-1));
measureddamageph =  log( pHsubset.Viabilty.^(-1));

damagepredictiontemp= exp(popttemp.logA ) * temperaturesubset.Time.*exp(-(popttemp.Ea0                                   )./(GasConst * temperaturesubset.Temperature));
measureddamagetemp =  log( temperaturesubset.Viabilty.^(-1));

% get Rsquared
mdl = fitlm(measureddamage ,damageprediction)
handlefit=figure(5)
plot(mdl)
xlabel( 'measured damage')
ylabel( 'predicted damage')
title(sprintf('R^2=%f, A=%9.2e, Ea0=%9.2e, Ea1=%9.2e',mdl.Rsquared.Ordinary,exp(popt.logA ) ,popt.Ea0,popt.Ea1))
saveas(handlefit,'ArrheniusFit','png')

% get Rsquared
mdlph = fitlm(measureddamageph ,damagepredictionph)
handlefit=figure(6)
plot(mdlph)
xlabel( 'measured damage')
ylabel( 'predicted damage')
title(sprintf('R^2=%f, A=%9.2e, Ea0=%9.2e, Ea1=%9.2e',mdlph.Rsquared.Ordinary,exp(poptph.logA ) ,poptph.Ea0,poptph.Ea1))
saveas(handlefit,'ArrheniusFitph','png')

% get Rsquared
mdltemp = fitlm(measureddamagetemp ,damagepredictiontemp)
handlefit=figure(7)
plot(mdltemp)
xlabel( 'measured damage')
ylabel( 'predicted damage')
title(sprintf('R^2=%f, A=%9.2e, Ea0=%9.2e',mdltemp.Rsquared.Ordinary,exp(popttemp.logA ) ,popttemp.Ea0))
saveas(handlefit,'ArrheniusFittemp','png')
