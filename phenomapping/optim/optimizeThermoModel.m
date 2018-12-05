function [sol] = optimizeThermoModel(tModel,minNorm,solver,time) %,UseIndConst,IndConstrNames)
% this function solves a TFBA problem using either Gurobi (through
% Gurobi_mex) or CPLEX
%

if (nargin < 2)
    minNorm = 0;
end
if (nargin < 3)
    solver = 'cplex'; %by default, another solver could be implemented here
end
if (nargin < 4)
    time = [];
end

sol = solveTFAmodelCplex(tModel,time);

if minNorm
    model4minNorm = tModel;
    model4minNorm.var_lb(model4minNorm.f==1) = sol.val;
    FRidx = getAllVar(model4minNorm,{'F','R'});
    f = zeros(numel(model4minNorm.f),1);
    f(FRidx) = 1;
    model4minNorm.f = f;
    model4minNorm.objtype = 1;
    sol = solveTFAmodelCplex(model4minNorm);
end

if ~isempty(sol.val) && ~isnan(sol.val)
    sol.stat = 1;
else
    sol.stat = 0;
    sol.val = NaN;
end