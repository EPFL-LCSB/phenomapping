function [pairsgenes, grRateEffectpairs] = extractDoubleGeneMatrix(model, grRateKOmatrix, geneList)
% Analyzes output of thermoDoubleGeneDeletion.m
%
% USAGE:
%
%    [pairsgenes, grRateEffectpairs] = extractDoubleGeneMatrix(model, grRateKOmatrix, genes)
%
% INPUT:
%    model:           TFA model structure including gene-reaction associations
%    grRateKOmatrix:  Matrix with the growth value upon double knockout of
%                     genes in geneList
% OPTIONAL INPUTS:
%    geneList:        Genes tested for doubleKO
%
% OUTPUTS:
%    pairsgenes:            Summary of gene pairs knocked out
%    grRateEffectpairs:     Growth value upon double knockout of gene
%                           paris in pairsgenes
%
% .. Author:
%       - Anush Chiappino-Pepe 31/8/17
%


if nargin < 3
    geneList = model.genes(:,1);
else
    [~,indx] = ismember(geneList, model.genes);
    if min(indx)==0
        error('not all genes were found');
    end
end

grRateKOmatrix = round(grRateKOmatrix,7);
rigthM = zeros(length(geneList),1);
leftM = zeros(length(geneList),1);

for i = 1:length(geneList)
    rigthM(i) = all((grRateKOmatrix(i,1:i) == grRateKOmatrix(1:i,i)'));
    leftM(i) = all((grRateKOmatrix(i,1:end+1-i) == grRateKOmatrix(1:end+1-i,end+1-i)'));
end

tmpR = tril(ones(length(geneList),length(geneList)));
tmpL = fliplr(tril(ones(length(geneList),length(geneList))));

if all(rigthM) || all((1000*tmpR(tmpR==1) == grRateKOmatrix(tmpR==1)))
    discard = tril(1000*ones(length(geneList),length(geneList)));
elseif all(leftM) || all((1000*tmpL(tmpL==1) == grRateKOmatrix(tmpL==1)))
    discard = fliplr(tril(1000*ones(length(geneList),length(geneList))));
else
    error('no simmetry found in the grRateKOmatrix')
end

% pairsgenes = cell(length(geneList)*length(geneList),2);
% grRateEffectpairs = zeros(length(geneList)*length(geneList),1);
pairsgenes = {};
grRateEffectpairs = [];

for i = 1:length(geneList)
%         pairsgenes(1+(i-1)*length(geneList):i*length(geneList),1) = geneList(i);
    for j = 1:length(geneList)
%                 pairsgenes(j+(i-1)*length(geneList),2) = geneList(j);
        pairsgenes(end+1,1) = geneList(i);
        pairsgenes(end,2) = geneList(j);
        if discard(i,j)==1000
            grRateEffectpairs(end+1,1) = 1000;
%             grRateEffectpairs(j+(i-1)*length(geneList),1) = 1000;
        else
            grRateEffectpairs(end+1,1) = grRateKOmatrix(i,j);
%             grRateEffectpairs(j+(i-1)*length(geneList),1) = grRateKOmatrix(i,j);
        end
    end
end