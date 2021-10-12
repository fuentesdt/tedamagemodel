function DamageDiff   = viability_model(Molarity, parnames, parvals)
global verbose;
% y = line_model(x, parnames, parvals)
%
% This function will return a line given by the equation y = mx + b, where
% m is the line's gradient and b is its y-intersept. The
% input parameters are:
%   x - the x values at which y will be calculated
%   parnames - a cell array containing the parameters names. These can be
%       in any order, but must include the following parameters:
%           {'m', 'b'}
%   parvals - a cell array containing the values of the parameters given in
%       parnames. These must be in the same order as in parnames.
%
%--------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
%--------------------------------------------------------------------------

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length!');
end

nparams = lpn;

% extract parameter values
for ii=1:nparams
  switch parnames{ii}
    case 'm'
      m = parvals{ii};
    case 'b'
      b = parvals{ii};
  end
end

% calculate y-values
%y = m*x + b;

% evaluate objective function
%solution =  [  6.28e5 100.] 
Ea = b +   Molarity   * m;
logA = 3.8e-4 * Ea - 9.36;
Deltat37degC = 24*60*60; % 24hr
Deltat43degC =  1*60*60; % 1hr
GasConst = 8.314;
%Damage37degC  = Deltat37degC *  exp(logA).*  exp(-Ea / GasConst / (37 +273)) ; 
%Damage43degC  = Deltat43degC *  exp(logA).*  exp(-Ea / GasConst / (43 +273)) + (Deltat37degC-Deltat43degC) * exp(logA).* exp(-Ea / GasConst / (37 +273)) ; 
DamageDiff  = Deltat43degC *  ( exp(logA-Ea / GasConst / (43 +273))  -  exp(logA-Ea / GasConst / (37 +273)) ) ; 
%Objective = sum((DamageDiff - logviabilityratio).^2)
    if verbose
        fprintf(1, 'b: %.5e, b = %.5e, DamageDiff  = %.5e %.5e %.5e %.5e %.5e \n', b,m,DamageDiff);
    end

return
