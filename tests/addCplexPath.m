function [solverOK,path_found] = addCplexPath(cplexPath)

if (nargin < 1)
    cplexPath = [];
end

while isempty(cplexPath)
    cplexPath = input('Please provide your cplex path and press enter\n... ','s');
end

addpath(genpath(cplexPath))
[path_found,~] = which('cplex');
solverOK = ~isempty(path_found);

if ~solverOK
    warning('you need to provide a valid path to the cplex solver')
end

% For now only cplex is available as a solver for FBA problems in cobra for
% PhenoMapping (see ext/cobratoolbox/solvers)
if ~changeCobraSolver('cplex_direct','LP')
    warning('you need to install cplex or add the solver functions from the cobratoolbox to the path')
end