%% Combined Thermal and Osmotic Stress Damage Model
% opts.format='pdf'; opts.outputDir='.'; publish('tcadamage.m',opts);
clear all
close all



% After 1h at 43degC the cells were returned to 37degC. After 24h of treatment initiation, cells were rinsed with PBS to
% removed the salt from the cultures and returned to standard culture conditions.

deltat = 10
time = [ 0:deltat:600]; % 5 mins
temperature = 37 + (80 - 37) * (1 - exp(- .002 * time));
figure(1)
plot(time,temperature)
xlabel('time')
ylabel('temperature')
% Arrhenius Model expect temperature in Kelvin
tempkelvin = temperature + 273;

% Henriques Arrhenius Parameters
HenriquesEa =  6.28e5
HenriquesA  = 3.1e98
GasConst = 8.314
log(HenriquesA ) -( 3.8e-4 * HenriquesEa - 9.36)

% Following Pearce Paper - 
%     Pearce, John A. "Comparative analysis of mathematical models
%     of cell death and thermal damage processes." International 
%     Journal of Hyperthermia  (2013): 262-280. 
% The Arrhenius damage parameters follow a log-linear relationship
Ea =  [4:2:16] * 1.e5;
logA = 3.8e-4 * Ea - 9.36;

% Plot expected damage for increaseing Activation Energy, Ea
% Notice that damage occurs faster with increasing activation energy
figure(2)
plotcol = [ 'b' ,'g' ,'r' ,'c' ,'m' ,'y' ,'k' ];
hold on
for iii = 1:length(Ea)
  damage = exp(logA(iii)) * deltat* exp(-Ea(iii) / GasConst *tempkelvin.^(-1)) ;
  damage = cumsum(damage);
  fractionalconversion = 1 - exp(- damage ) ;
  plot(time,fractionalconversion, plotcol(iii) )
end
xlabel('time')
ylabel('FC')

% consider a second temperature that max out at 47degC
temperature = min(37 + (80 - 37) * (1 - exp(- .002 * time)),47);
tempkelvin = temperature + 273;
figure(3)
plot(time,temperature)
xlabel('time')
ylabel('temperature')
% at the same time, we will consider a steadly increasing salt concentration
saltconcentration = [0:.2:1]
% assume activation energy a function of the salt concentration
Ea = 4e5 +  saltconcentration * (16e5-4e5);
logA = 3.8e-4 * Ea - 9.36;

% Plot expected damage for increasing Activation Energy, Ea
% Notice that damage occurs faster with increasing activation energy/salt concentration
% At these lower temperature the damage still increases with increasing salt presence
figure(4)
plotcol = [ 'b' ,'g' ,'r' ,'c' ,'m' ,'y' ,'k' ];
hold on
for iii = 1:length(Ea)
  damage = exp(logA(iii)) * deltat* exp(-Ea(iii) / GasConst *tempkelvin.^(-1)) ;
  damage = cumsum(damage);
  fractionalconversion = 1 - exp(- damage ) ;
  plot(time,fractionalconversion, plotcol(iii) )
end
xlabel('time')
ylabel('FC')
legend('c=0.0','c=0.2','c=0.4','c=0.6','c=0.8','c=1.0')


% consider another temperature that max out at 42degC
temperature = min(37 + (80 - 37) * (1 - exp(- .002 * time)),42);
tempkelvin = temperature + 273;
figure(5)
plot(time,temperature)
xlabel('time')
ylabel('temperature')
% at the same time, we will consider a steadly increasing salt concentration
saltconcentration = [0:.2:1]
% assume activation energy a function of the salt concentration
Ea = 4e5 +  saltconcentration * (16e5-4e5);
logA = 3.8e-4 * Ea - 9.36;

% Plot expected damage for increasing Activation Energy, Ea
% Notice that damage occurs faster with increasing activation energy/salt concentration
% At these lower temperature the damage still increases with increasing salt presence
figure(6)
plotcol = [ 'b' ,'g' ,'r' ,'c' ,'m' ,'y' ,'k' ];
hold on
for iii = 1:length(Ea)
  damage = exp(logA(iii)) * deltat* exp(-Ea(iii) / GasConst *tempkelvin.^(-1)) ;
  damage = cumsum(damage);
  fractionalconversion = 1 - exp(- damage ) ;
  plot(time,fractionalconversion, plotcol(iii) )
end
xlabel('time')
ylabel('FC')
ylim([0 1])
legend('c=0.0','c=0.2','c=0.4','c=0.6','c=0.8','c=1.0')
