% addpath /rsrch1/ip/dtfuentes/github/matlabmultinest/general
% addpath /rsrch1/ip/dtfuentes/github/matlabmultinest/src
% example: (from hogg et al., 1008.4686)
% fit a line to data omitting outlier points
clear all
close all

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

% read sample data for fitting a line
%data = readdata_line;

% omit data corresponding to outliers (first 4 points)
%data{1} = data{1}(5:end); % x_i
%data{2} = data{2}(5:end); % y_i
%data{3} = data{3}(5:end); % sigma_yi

% given data
Viability37degC = .01*[100	101	95	87	22];
Viability43degC = .01*[84	84	82	66	12];
FC37degC = 1- Viability37degC;
FC43degC = 1- Viability43degC;
Molarity        =    [  0	50	100	200	400]; % mM
logviabilityratio =  - log((Viability37degC.^-1).*Viability43degC )

data{1} = Molarity'      ; % x_i
data{2} = logviabilityratio' ; % y_i
% convert sigmas to a covariance matrix and reassign to data{3}
C = diag(ones(length(Molarity),1))
data{3} = C;

% define nested sampling parameters
Nlive = 500;
%Nmcmc = input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
Nmcmc = 0
tolerance = .1;
%tolerance = 1.e10;
likelihood = @logL_gaussian;
model = @viability_model;
prior = {'m', 'uniform', 0, 10000, 'fixed'; ...
         'b', 'uniform', 6.28e1 , 6.28e6 , 'fixed'};
extraparams = {}; % no extra parameters beyond what's in the prior

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

% plot posterior distributions
wp = [1];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [2];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [1 2];
posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});

