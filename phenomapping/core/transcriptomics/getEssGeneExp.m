function [addEssGenesExp, essGenesExp, grRateExp] = getEssGeneExp...
    (expModel, sol, model, minObj, essTFAref)
% Identify essential genes per gene expression profile in a model that
% integrates expression constraints from TEX-FBA
%
% USAGE:
%
%       [grRateExp, essGenesExp, addEssGenesExp] = getEssGeneExp(expModel, sol, model, minObj, essTFAref)
%
% INPUTS:
%    expModel:        Model with gene expression constraints from TEX-FBA
%                     (obtained with function: getExpProfile.m)
%    sol:             sol obtained in getExpProfile.m from TEX-FBA
%                     where alternative sols are saved
%    model:           TFA model structure (the "model" input to 
%                     integrateGeneExp.m)
%
% OPTIONAL INPUTS:
%    minObj:          Minimal objective value desired (default = 0.05)
%    essTFAref:       Essentiality without gene expression (default = 
%                     tested with thermoSingleGeneDeletion.m)
%
% OUTPUTS:
%    addEssGenesExp:  Essential genes with gene expression on the top of 
%                     essential genes without gene expression
%    essGenesExp:     Essential genes per expression profile
%    grRateExp:       Predicted growth upon knockout of genes per
%                     expression profile
%
% .. Author:
% Anush Chiappino-Pepe 2017
% 

if (nargin < 4) || isempty(minObj)
    minObj = 0.05;
end
if (nargin < 5)
    fprintf('TFA essentiality will be performed\n');
    [~, grRateKOGeneTFA] = thermoSingleGeneDeletion(model, 'TFA', ...
        model.genes, 0, 0, 0, 0, model.indNF);
    grRateKOGeneTFA(isnan(grRateKOGeneTFA)) = 0;
    essTFAref = model.genes(grRateKOGeneTFA < minObj);
end

indUPDOWN = [expModel.ExpInfo.indUP; expModel.ExpInfo.indDOWN];

essGenesExp = cell(size(sol.z_matrix,2),1);
addEssGenesExp = cell(size(sol.z_matrix,2),1);
grRateExp = cell(size(sol.z_matrix,2),1);

posaltmcs = find(sol.store_obj>max(sol.store_obj)-0.5);
 
for i = 1:length(posaltmcs)
    expModelAlt = integrateConsProfile(expModel, sol, indUPDOWN, ...
        posaltmcs(i));
    [~, grRateExp{i}] = thermoSingleGeneDeletion(expModelAlt,...
        'TFA', expModelAlt.genes, 0, 0, 0, 0, expModelAlt.indNF);
    grRateExp{i}(isnan(grRateExp{i})) = 0;
    essGenesExp{i} = expModel.genes(grRateExp{i} < minObj);
    addEssGenesExp{i} = essGenesExp{i}(~ismember(essGenesExp{i}, ...
        essTFAref));
end