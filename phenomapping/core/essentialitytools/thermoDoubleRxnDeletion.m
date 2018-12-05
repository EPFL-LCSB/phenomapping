function [grRateKOmatrix,grRatioDble,grRateKO,grRateWT] = thermoDoubleRxnDeletion(model, rxnList1, rxnList2)
% Performs double rxn deletion analysis using TFA
%
% USAGE:
%
%    [grRateKOmatrix] = thermoDoublerxnDeletion(model, rxnList1, rxnList2)
%
% INPUT:
%    model:           TFA model structure including rxn-reaction associations
%
% OPTIONAL INPUTS:
%    rxnList:        List of rxns to be deleted (default = all rxns)
%    essThr:          Essentiality threshold for single rxn deletion
%
% OUTPUTS:
%    grRateKOmatrix:  Matrix with the growth value upon double knockout of
%                     rxns in rxnList
%    rxnList:        rxns tested for doubleKO
%
% .. Author:
%       - Anush Chiappino-Pepe 31/8/17
%
% Note: analyze the grRateKOmatrix with extractDoublerxnMatrix.m

differentSetsFlag = false;

if (nargin < 2)
    rxnList1 = model.rxns;
else
    if (isempty(rxnList1))
        rxnList1 = model.rxns;
    end
end
if (nargin < 3)
    rxnList2 = rxnList1;
else
    if (isempty(rxnList2))
        rxnList2 = rxnList1;
    else
        differentSetsFlag = true;
    end
end

indNF = getAllVar(model,{'NF'});
if isempty(indNF)
    model = addNetFluxVariablesNEW(model);
end

% if not provided, calculate essential rxns with tfba (single KO)
% Extract the rxns that are not essential with single thermo rxn knockout
fprintf('Single deletion analysis to remove lethal rxns\n');
[singleRatio1,singleRate1,grRateWT] = thermoSingleRxnDeletion(model,'TFA',rxnList1);
singleLethal1 = (singleRatio1 < 1e-9);
rxnListOrig1 = rxnList1;
rxnListOrig2 = rxnList1;
rxnList1 = rxnList1(~singleLethal1);
singleRate = singleRate1(~singleLethal1);
[~,listMap1] = ismember(rxnListOrig1,rxnList1);
fprintf('%d non-lethal rxns\n',length(rxnList1));

% Repeat the analysis for the second set of rxns
if (differentSetsFlag)
    fprintf('Single deletion analysis to remove lethal rxns from rxn set 2\n');
    [singleRatio2] = thermoSingleRxnDeletion(model,'TFA',rxnList2);
    singleLethal2 = (singleRatio2 < 1e-9);
    rxnListOrig2 = rxnList2;
    rxnList2 = rxnList2(~singleLethal2);
    [~,listMap2] = ismember(rxnListOrig2,rxnList2);
    fprintf('%d non-lethal rxns\n',length(rxnList2));
else
    rxnList2 = rxnList1;
    listMap2 = listMap1;
end

nDelrxns1 = length(rxnList1);
nDelrxns2 = length(rxnList2);

grRateKO = ones(nDelrxns1,nDelrxns2)*grRateWT;
allGrRateKO = 100*ones(nDelrxns1,nDelrxns2);

if (differentSetsFlag)
    nTotalPairs = nDelrxns1*nDelrxns2;
else
    nTotalPairs = nDelrxns1*(nDelrxns1-1)/2;
end

% grRateKOmatrix will be a simetric matrix; avoid calculating twice the
% same by discarding elements = 1000 (lower triangular matrix)
grRateKOmatrix = tril(1000*ones(length(rxnList1),length(rxnList1)));
grRateWT = optimizeThermoModel(model);
grRateWT = grRateWT.val;

% For each of these rxns create a model without it and perform single rxn
% deletion
% Study the combined effect

delCounter = 0;
fprintf('Double rxn deletion analysis\n');
fprintf('Total of %d pairs to analyze\n',nTotalPairs);
showprogress(0,'Double rxn deletion analysis in progress ...');
t = cputime;
fprintf('Perc complete\tCPU time\n');
for rxnNo1 = 1:nDelrxns1
    [~ ,rxnID1] = ismember(rxnList1{rxnNo1},model.rxns);
    if (~differentSetsFlag)
        grRateKO(rxnNo1,rxnNo1) = singleRate(rxnNo1);
        initID = rxnNo1+1;
    else
        initID = 1;
    end
    for rxnNo2 = initID:nDelrxns2
        delCounter = delCounter + 1;
        if (mod(delCounter,10) == 0)
            showprogress(delCounter/nTotalPairs);
        end
        if (mod(delCounter,100) == 0)
            fprintf('%5.2f\t%8.1f\n',100*delCounter/nTotalPairs,cputime-t);
        end
        % Save results every 1000 steps
        if (mod(delCounter,1000) == 0)
            save doublerxnDeletionTmp.mat grRateKO grRateKOmatrix
        end
        [~ ,rxnID2] = ismember(rxnList2{rxnNo2},model.rxns);
        % Use FBA/MOMA/lMOMA to calculate deletion strain growth rate
        constrRxnInd = [rxnID1; rxnID2];
        [~ , ba] = ismember(strcat('NF_', model.rxns(constrRxnInd)), model.varNames);
        dmodel = model;
        dmodel.var_ub(ba) = 0;
        dmodel.var_lb(ba) = 0;
        solKO = optimizeThermoModel(dmodel);
        
        if (solKO.stat == 1)
            grRateKO(rxnNo1,rxnNo2) = solKO.val;
            grRateKO(rxnNo2,rxnNo1) = solKO.val;
            grRateKOmatrix(rxnNo1,rxnNo2) = solKO.val;
        else
            grRateKO(rxnNo1,rxnNo2) = NaN;
            grRateKO(rxnNo2,rxnNo1) = NaN;
            grRateKOmatrix(rxnNo1,rxnNo2) = NaN;
        end
        if (differentSetsFlag)
            grRateKO(rxnNo2,rxnNo1) = grRateKO(rxnNo1,rxnNo2);
        end
    end
end

% Reconstruct the entire matrix
for i = 1:length(rxnListOrig1)
    for j = 1:length(rxnListOrig2)
        if (listMap1(i) > 0 && listMap2(j) > 0)
            allGrRateKO(i,j) = grRateKO(listMap1(i),listMap2(j));
        else
            allGrRateKO(i,j) = 0;
        end
    end
end

grRatioDble = allGrRateKO/grRateWT;

grRateKO = allGrRateKO;
