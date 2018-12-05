function [blockedRxns] = essenEvalRxns(model,rxns,rxnGPRs,genes,geneRatio)
% Identify reactions to delete if we assume that the lowly expressed genes
% are "knocked out"
%
% [blockedRxns] = essenEvalRxns(model,rxns,rxnGPRs,genes,geneRatio)
%
%INPUT
% model             TFA model structure
%
%OPTIONAL INPUTS
% rxns              Include objective rxns in the exchange rxn set (1) or not (0)
%                   (Default = false)
% irrevFlag         Model is in irreversible format (1) or not
%                   (Default = false)
%
%OUTPUTS
% blockedRxns       Reactions to delete (associated to lowly expressed
%                   genes)
%                   model is exchange or not 
% selUpt            Boolean vector indicating whether each reaction in
%                   model is nutrient uptake or not
%
% Exchange reactions only have one non-zero (+1/-1) element in the 
% corresponding column of the stoichiometric matrix. Uptake reactions are 
% exchange reactions are exchange reactions with negative lower bounds.
%
% 10/14/05 Markus Herrgard

% MUmodel is TFA fomat model
% geneIdx: gene indexes in the model
% geneRatio: coreesponding gene ratios
% output

geneRatio_map = containers.Map(genes,geneRatio);
rxnExp = nan(numel(rxns),1);
for r = 1:numel(rxns)
    rxnExp(r) = evalEachRxn(model,rxns{r},rxnGPRs{r},geneRatio_map,...
        @min,@max);
end
blockedRxns = rxns(rxnExp==0);
end