function [blockedRxns] = essenEvalRxns(model, rxns, rxnGPRs, genes,...
    geneRatio)
% Identify reactions associated to lowly expressed genes
%
% [blockedRxns] = essenEvalRxns(model, rxns, rxnGPRs, genes, geneRatio)
%
% INPUT
% model             TFA model structure
% rxns              Include objective rxns in the exchange rxn set (1) or not (0)
%                   (Default = false)
% rxnGPRs           Model is in irreversible format (1) or not
%                   (Default = false)
%
% OUTPUTS
% blockedRxns       Reactions associated to lowly expressed genes
% 
%
% 2016 Vikash Pandey

geneRatio_map = containers.Map(genes,geneRatio);
rxnExp = nan(numel(rxns),1);
for r = 1:numel(rxns)
    rxnExp(r) = evalEachRxn(model,rxns{r},rxnGPRs{r},geneRatio_map,...
        @min,@max);
end
blockedRxns = rxns(rxnExp==0);
end