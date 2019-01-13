function botRxnLevels = linkEssGene2Exp(expModel, essExp, sol, indInt, ...
    minObj, NumAlt, indNF, filename)
% Link essential gene to underlying reaction levels
%
% botRxnLevels = linkEssGene2Exp(expModel, essExp, sol, indInt, minObj, NumAlt, indNF, filename)
%
%
% INPUT
% expModel          TEX-FBA expModel structure with objective (in 
%                   expModel.f) set for the objective variable of study
% essExp:           Essentiality per expression profile as output of 
%                   getEssGeneExp.m 
% sol:              sol obtained in getExpProfile.m from TEX-FBA
%                   where alternative sols are saved
%
% OPTIONAL INPUTS:
% indInt            Indexes of integer variables for analysis, in TEX-FBA  
%                   these are the UP and DOWN rxns (default = all of them)
% minObj:           Minimal objective value desired (default = 0.05)
% NumAlt:           Maximum number of alternatives to get (default = 1)
% indNF:            Indexes of net fluxes of all reactions (default = 
%                   get them here)
% filename:         Name used to save DPs (default =
%                   'PhenoMappingTranscriptomics')
%
% OUTPUTS
% botRxnLevels      Bottleneck rxn levels in TEX-FBA model
%
%
% .. Author:
% Anush Chiappino-Pepe 2017
%

if (nargin < 4) || isempty(indInt)
    indInt = [expModel.ExpInfo.indUP; expModel.ExpInfo.indDOWN];
end
if (nargin < 5) || isempty(minObj)
    minObj = 0.05;
end
if (nargin < 6) || isempty(NumAlt)
    NumAlt = 1;
end
if (nargin < 7) || isempty(NumAlt)
    indNF = getAllVar(expModel, {'NF'});
end
if (nargin < 8) || isempty(filename)
    filename = 'PhenoMappingTranscriptomics';
end

posaltmcs = find(sol.store_obj>max(sol.store_obj)-0.5);
botRxnLevels = cell(length(posaltmcs),1);

for sAi = 1:length(posaltmcs) % select alternative to integrate   
    indInti = indInt(sol.sol_matrix(indInt,posaltmcs(sAi))>0.9);
    for i = 1:length(essExp{sAi})
        modeli = thermoDeleteModelGenes(expModel, essExp{sAi}(i));
        rxnsToRelax = getBotNeckRxns(modeli, indInti, minObj, NumAlt,...
            indNF, filename);
        botRxnLevels{sAi}(i,:) = rxnsToRelax;
    end
end