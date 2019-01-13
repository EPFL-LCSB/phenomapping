function [minmax, description, report] = checkDirectionality(model, ...
    minmax, tagFBA, perObj, precision)
% Get reaction directionalities for the FVA output
%
% USAGE:
%
%    [description, report] = checkDirectionality(model, minmax, tagFBA, perObj)
%
% INPUT:
%    model:           FBA/TFA model structure
%
% OPTIONAL INPUTS:
%    minmax:          FVA result with lb in 1st column and ub in 2nd column
%                     (default = FVA for all rxns in the model)
%    tagFBA:          True for FBA, false for TFA (default = true)
%    perObj:          % of optimal objective value to fix for a FVA
%                     (default = 90)
%    precision:       Rounding minmax output (default = 9 decimals)
%
% OUTPUTS:
%    minmax:          Minmax used for directionality analysis
%    description:     Directionality classification for the minmax given
%    report:          Nb of reactions of each class in the minmax given
%
% .. Author:
% Anush Chiappino-Pepe 2014
%

if (nargin < 2)
    minmax = [];
end
if (nargin < 3) || isempty(tagFBA)
    tagFBA = 1;
end
if (nargin < 4) || isempty(perObj)
    perObj = 0.9;
end
if (nargin < 5) || isempty(precision)
    precision = 9;
end

% require some growth
if tagFBA
    sol = solveFBAmodelCplex(model);
    model.lb(model.c==1) = sol.f*perObj;
else
    sol = optimizeThermoModel(model);
    model.var_lb(model.f==1) = sol.val*perObj;
end

% get minmax if not there and round
if isempty(minmax) && tagFBA
    minmax = runMinMax(model);
elseif isempty(minmax) && ~tagFBA
    minmax = runTMinMax(model,strcat('NF_',model.rxns),1);
end
minmax = round(minmax, precision);

% reaction directionality based on FVA/minmax
description = cell(size(minmax,1),1);
for i = 1:size(minmax,1)
    % negative lb: it can only be bidirectional or backwards
    if minmax(i,1) < 0 
        if minmax(i,2) > 0 % positive ub
            description(i) = {'bidirectional'};
        elseif minmax(i,2) <= 0 % negative ub
            description(i) = {'backwards'};
        end
    % positive or zero ub: it can only be forwards or blocked
    elseif minmax(i,1) >= 0
        if minmax(i,2) > 0 % positive ub
            description(i) = {'forwards'};
        elseif minmax(i,2) <= 0 % zero ub
            description(i) = {'blocked'};
        end
    end
end

% report of reactions directionalities
report = {'bidirectional';'forwards';'backwards';'blocked'};
report{1,2} = num2str(sum(ismember(description,'bidirectional')));
report{2,2} = num2str(sum(ismember(description,'forwards')));
report{3,2} = num2str(sum(ismember(description,'backwards')));
report{4,2} = num2str(sum(ismember(description,'blocked')));
