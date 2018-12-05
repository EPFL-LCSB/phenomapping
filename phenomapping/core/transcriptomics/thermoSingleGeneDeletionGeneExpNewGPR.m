function [objvals, essgenes, allTasksImpact, genesToRelax, model1] = thermoSingleGeneDeletionGeneExpNewGPR...
(model, grRate, essThr, flagTasks, geneList, flagBotNeck, indObj,numAlt,indNF)
% Determine essential genes for a model with a gene expression profile
% integrated
%
% USAGE:
%
%       [objvals, essgenes, allTasksImpact] = thermoSingleGeneDeletionGeneExpNewGPR(model, grRate, essThr, flagTasks)
%
% INPUTS:
%    model:           model with TFA structure and a gene expression
%                     profile integrated (from integrateConsProfile.m)
%
%
% OPTIONAL INPUTS:
%    grRate:          growth rate desired (default = solution from 
%                     optimization of current model)
%    essThr:          essenthiality threshold (default = 0.1)
%    flagTasks:       perform analysis of essential tasks (default = 0)
%
% OUTPUTS:
%    objvals:         grRate upong gene deletion
%    essgenes:        essential genes identified
%    allTasksImpact:  metabolic tasks impacted upon knockout of essential gene
%
% Anush Chiappino-Pepe 2017
%
if (nargin < 2)
    solnew = optimizeThermoModel(model);
    if isempty(solnew.val) || isnan(solnew.val)
        error('model is not feasible')
    end
    grRate = solnew.val;
end
if (nargin < 3)
    essThr = 0.1;
end
if (nargin < 4)
    flagTasks = 0;
end
if (nargin < 5)
    geneList = model.genes;
end
if (nargin < 6)
    flagBotNeck = 0;
end
if (nargin < 7)
    indObj = [model.ExpInfo.indUP; model.ExpInfo.indDOWN];
end
if (nargin < 8)
    numAlt = 1;
end
if (nargin < 9)
    indNF = [];
end

genes = model.genesnewgrRules;
genes = genes(ismember(genes,geneList));
genesToRelax = cell(length(genes),1);

rxns = model.rxnsnewgrRules;
newgrRules = model.newgrRules;

impactTasks = cell(numel(genes),3);
alwaysOne = ones(numel(genes),1);
objvals = ones(numel(genes),1)*1000;

mingrRate = essThr*grRate;
for i = 1:numel(genes)
    geneRatio = alwaysOne;
    geneRatio(i) = 0;
    [blockedRxns] = essenEvalRxns(model,rxns,newgrRules,genes,geneRatio);
    [model1] = thermoDeleteModelRxns(model,blockedRxns);
    sol = optimizeThermoModel(model1);
    if ~isempty(sol.x)
        objvals(i) = sol.val;
        if (sol.val < mingrRate) && flagTasks
            [impactTasks{i,1}, impactTasks{i,2}, impactTasks{i,3}] = checkBBBTasks(model1,'TFA',essThr,[],0,0,grRate);
        elseif (sol.val < mingrRate) && flagBotNeck
            [~,genesToRelax{i,1}] = getBotGenes(model1,indObj,mingrRate,numAlt,indNF);
        end
    end
end
% objvals2 = objvals;
objvals(objvals==1000) = 0;
essgenes = genes(objvals<mingrRate)';

if flagTasks
    allTasksImpact = extractTasksImpact(impactTasks(ismember(genes,essgenes),3), essgenes);
else
    allTasksImpact = {};
end

