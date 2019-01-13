function [rxnsToRelax, rxnsToKeep] = getBotNeckRxns(expModel, indInt, ...
    minObj, NumAlt, indNF, filename)
% Identify reaction levels to relax and make model feasible
%
% [rxnsToRelax, rxnsToKeep] = getBotNeckRxns(expModel, indInt, minObj, NumAlt, indNF, filename)
%
%
% INPUT
% expModel          TEX-FBA expModel structure with objective (in 
%                   expModel.f) set for the objective variable of study
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
% rxnsToRelax       Rxn levels that cannot be kept in a feasible model or
%                   bottleneck rxn levels
% rxnsToKeep        Rxn levels that can be kept in a feasible model
%
%
% .. Author:
% Anush Chiappino-Pepe 2017
% 

if (nargin < 2) || isempty(indInt)
    indInt = [expModel.ExpInfo.indUP; expModel.ExpInfo.indDOWN];
end
if (nargin < 3) || isempty(minObj)
    minObj = 0.05;
end
if (nargin < 4) || isempty(NumAlt)
    NumAlt = 1;
end
if (nargin < 5) || isempty(NumAlt)
    indNF = getAllVar(expModel, {'NF'});
end
if (nargin < 6) || isempty(filename)
    filename = 'PhenoMappingTranscriptomics_DPs';
end

rxnsToKeep = cell(NumAlt,1);
rxnsToRelax = cell(NumAlt,1);

% remove maximum consistency score (MCS) requirement
expModel.var_lb(ismember(expModel.varNames,{'Total_UpsAndDowns'})) = 0;

% set growth requirement
expModel.var_lb(expModel.f==1) = minObj;

% set objective in rxns that caused essentiality
expModel.f = zeros(length(expModel.varNames),1);
expModel.f(indInt) = 1;
expModel.objtype = -1; %max

solution = getExpProfile(expModel, NumAlt, indInt, indNF, filename);

rxns = expModel.varNames(indInt);
if ~isempty(solution.sol_matrix)
    DP = solution.sol_matrix(indInt,:);
else
    DP =[];
    rxnsToRelax = [];
    rxnsToKeep = rxns;
end

for i = 1:size(DP,2)
    rxnsToKeep{i} = rxns(DP(:,i)>0.9);
    rxnsToRelax{i} = rxns(DP(:,i)<0.1);
end
if size(DP,2)<NumAlt
    rxnsToKeep = rxnsToKeep(1:size(DP,2));
    rxnsToRelax = rxnsToRelax(1:size(DP,2));
end
