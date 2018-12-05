function [pairsrxns, grRateEffectpairs] = extractDoubleRxnMatrix(model, grRateKOmatrix, rxnList)
% Analyzes output of thermoDoublerxnDeletion.m
%
% USAGE:
%
%    [pairsrxns, grRateEffectpairs] = extractDoublerxnMatrix(model, grRateKOmatrix, rxns)
%
% INPUT:
%    model:           TFA model structure including rxn-reaction associations
%    grRateKOmatrix:  Matrix with the growth value upon double knockout of
%                     rxns in rxnList
% OPTIONAL INPUTS:
%    rxnList:        rxns tested for doubleKO
%
% OUTPUTS:
%    pairsrxns:            Summary of rxn pairs knocked out
%    grRateEffectpairs:     Growth value upon double knockout of rxn
%                           paris in pairsrxns
%
% .. Author:
%       - Anush Chiappino-Pepe 31/8/17
%


if nargin < 3
    rxnList = model.rxns(:,1);
else
    [~,indx] = ismember(rxnList, model.rxns);
    if min(indx)==0
        error('not all rxns were found');
    end
end

grRateKOmatrix = round(grRateKOmatrix,7);
rigthM = zeros(length(rxnList),1);
leftM = zeros(length(rxnList),1);

for i = 1:length(rxnList)
    rigthM(i) = all((grRateKOmatrix(i,1:i) == grRateKOmatrix(1:i,i)'));
    leftM(i) = all((grRateKOmatrix(i,1:end+1-i) == grRateKOmatrix(1:end+1-i,end+1-i)'));
end

tmpR = tril(ones(length(rxnList),length(rxnList)));
tmpL = fliplr(tril(ones(length(rxnList),length(rxnList))));

if all(rigthM) || all((1000*tmpR(tmpR==1) == grRateKOmatrix(tmpR==1)))
    discard = tril(1000*ones(length(rxnList),length(rxnList)));
elseif all(leftM) || all((1000*tmpL(tmpL==1) == grRateKOmatrix(tmpL==1)))
    discard = fliplr(tril(1000*ones(length(rxnList),length(rxnList))));
else
    error('no simmetry found in the grRateKOmatrix')
end

% pairsrxns = cell(length(rxnList)*length(rxnList),2);
% grRateEffectpairs = zeros(length(rxnList)*length(rxnList),1);
pairsrxns = {};
grRateEffectpairs = [];

for i = 1:length(rxnList)
%         pairsrxns(1+(i-1)*length(rxnList):i*length(rxnList),1) = rxnList(i);
    for j = 1:length(rxnList)
%                 pairsrxns(j+(i-1)*length(rxnList),2) = rxnList(j);
        pairsrxns(end+1,1) = rxnList(i);
        pairsrxns(end,2) = rxnList(j);
        if discard(i,j)==1000
            grRateEffectpairs(end+1,1) = 1000;
%             grRateEffectpairs(j+(i-1)*length(rxnList),1) = 1000;
        else
            grRateEffectpairs(end+1,1) = grRateKOmatrix(i,j);
%             grRateEffectpairs(j+(i-1)*length(rxnList),1) = grRateKOmatrix(i,j);
        end
    end
end