function [bottleneckMets, modelnew, bottleneckMetNames] = getBotNeckMets ...
    (model, LCcons, grRate, NumAlt, time, tagMax, filename)
% Identify concentrations to relax to make model feasible
%
% USAGE:
%
%       [bottleneckMets, modelnew] = getBotNeckMets(model, LCcons, grRate, NumAlt)
%
% INPUTS:
%    model:           FEASIBLE model with TFA structure
%    LCcons:          Default structure to integrate concentrations in TFA
%
% OPTIONAL INPUTS:
%    grRate:          min growth to achieve (default = 0.007)
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%    tagMax:          True to generate alternative solution of minimal size
%                     only. False to generate up to the number of
%                     alternatives provided in NumAlt (default = false)
%    filename:        Name used to save DPs (default = 'PhenoMappingDPMin')
%
% OUTPUTS:
%    bottleneckMets:  Cell array with metabolites found in each alternative
%    modelnew:        Model with MILP formulation and integer cuts integrated
%
% .. Author:
% Meric Ataman & Anush Chiappino-Pepe 2015
% 

if (nargin < 3)
    grRate = 0.1*0.07;
end
if (nargin < 4)
    NumAlt = 100;
end
if (nargin < 5)
    time = [];
end
if (nargin < 6)
    tagMax = 0;
end
if (nargin < 7)
    filename = 'PhenoMappingDPMin';
end

conc = LCcons(:,1);
conc = conc(ismember(conc, model.varNames));
concBounds = LCcons(ismember(conc, model.varNames),2:3);
if iscell(concBounds)
    concBounds = cell2mat(concBounds);
end

indLCUSE = zeros(length(conc),1);
model.var_lb(model.f==1) = grRate;
model.f = zeros(length(model.f),1);
[~, lcrow] = ismember(conc, model.varNames); %LC constrained

fprintf('defining MILP problem\n');
intTag = {'LCUSE'}; % integer tag used in this MILP
[~,numVars] = size(model.A);
for i=1:length(conc)
    model.varNames(numVars+i,1) = strcat(intTag, '_', model.varNames(lcrow(i)));
    model.var_ub(numVars+i,1) = 1;
    model.var_lb(numVars+i,1) = 0;
    model.vartypes{numVars+i,1} = 'B';
    model.f(numVars+i,1)=1;
    indLCUSE(i,1)=numVars+i;
end

% define contraints for MILP
% LC + (ubLC + maxLC)*LCUSE < maxLC
% LCUSE=1 -> LC < ubLC
% LCUSE=0 -> LC < maxLC
[numCons,~] = size(model.A);
for i=1:length(conc)
    model.constraintNames{numCons+i} = ...
        strcat('LCUSE_relax_', model.varNames{lcrow(i)});
    model.rhs(numCons+i) = roundn(concBounds(i,2),-5);
    model.A(numCons+i, lcrow(i)) = 1;
    model.A(numCons+i, indLCUSE(i)) = ...
        roundn((-model.var_ub(lcrow(i))) + concBounds(i,2),-5);
    model.constraintType{numCons+i} = '<';
end

% LC + (lbLC + minLC)*LCUSE > minLC
% LCUSE=1 -> LC > lbLC
% LCUSE=0 -> LC > minLC
[numCons,~] = size(model.A);
for i=1:length(conc)
    model.constraintNames{numCons+i} = ...
        strcat('LCUSE_relax_2', model.varNames{lcrow(i)});
    model.rhs(numCons+i) = roundn(concBounds(i,1),-5);
    model.A(numCons+i, lcrow(i)) = 1;
    model.A(numCons+i, indLCUSE(i)) = ...
        roundn((-model.var_lb(lcrow(i))) + concBounds(i,1),-5);
    model.constraintType{numCons+i} = '>';
end

% minimization: try to implement all maxLC and minLC possible
fprintf('getting alternative solutions\n');
model.objtype = 1;
bottleneckMets = cell(1,length(NumAlt));
bottleneckMetNames = cell(1,length(NumAlt));
[cDPs, modelnew] = findDPMin(model, NumAlt, indLCUSE, time, tagMax, filename); % get alternative solutions
if isempty(cDPs)
    disp('no solution found');
else
    for i = 1:size(cDPs,2) % get metNames
        rxnsBFUSE = model.varNames(indLCUSE(cDPs(indLCUSE,i)>0.98));
        if isempty(rxnsBFUSE)
            error('no solution was found')
        else
            bottleneckMets{1,i} = strrep(rxnsBFUSE, strcat(intTag,'_LC_'), '');
            [~, tmp] = ismember(bottleneckMets{1,i}, model.mets);
            bottleneckMetNames{1,i} = model.metNames(tmp);
        end
    end
end