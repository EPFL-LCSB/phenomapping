function [genesToKeep,genesToRelax] = getBotGenes(model,indObj,grRate,numAlt,indNF,pathSave)

if (nargin < 2)
    indObj = [model.ExpInfo.indUP; model.ExpInfo.indDOWN];
end
if (nargin < 3)
    sol = optimizeThermoModel(model);
    grRate = sol.val;
end
if (nargin < 4)
    numAlt = 1;
end
if (nargin < 5)
    indNF = [];
end
if (nargin < 6)
    pathSave = '/Documents';
end

if isempty(indNF)
    indNF = getAllVar(model,{'NF'});
end

genesToKeep = cell(numAlt,1);
genesToRelax = cell(numAlt,1);

% remove CS requirement
model.var_lb(ismember(model.varNames,{'Total_UpsAndDowns'})) = 0;

% set growth requirement
model.var_lb(model.f==1) = grRate;

% set objective in rxns that caused essentiality
model.f = zeros(length(model.varNames),1);
model.f(indObj) = 1;
model.objtype = -1; %max

solution = getExpProfile(model, numAlt, indObj, indNF, pathSave);

DP = solution.sol_matrix(indObj,:);
rxns = model.varNames(indObj);

for i = 1:size(DP,2)
    genesToKeep{i} = rxns(DP(:,i)>0.9);
    genesToRelax{i} = rxns(DP(:,i)<0.1);
end
if size(DP,2)<numAlt
    genesToKeep = genesToKeep(1:size(DP,2));
    genesToRelax = genesToRelax(1:size(DP,2));
end
