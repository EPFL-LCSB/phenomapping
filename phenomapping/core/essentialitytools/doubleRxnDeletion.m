function [grRatioDble,grRateKO,grRateWT] = doubleRxnDeletion(model,method,rxnList1,rxnList2,verbFlag)
%doublerxnDeletion Performs double rxn deletion analysis using FBA, MOMA,
%or linear MOMA
%
% [grRatioDble,grRateKO,grRateWT] =
%     doubleRxnDeletion(model,method,rxnList1,rxnList2,verbFlag)
%
%INPUT
% model         COBRA model structure
%
%OPTIONAL INPUTS
% method        Either 'FBA' (default) 'MOMA', or 'lMOMA'
% rxnList1      List of query rxns to be deleted (default = all rxns)
% rxnList2      List of target rxns to be deleted (default = rxnList1)
% verbFlag      Verbose output (default = false)
%
%OUTPUTS
% grRatioDble   Computed growth rate ratio between double deletion strain 
%               and wild type
% grRateKO      Double deletion strain growth rates (1/h)
% grRateWT      Wild type growth rate (1/h)
%
% Markus Herrgard 8/8/06

if (nargin < 2)
    method = 'FBA';
end
if (nargin < 3)
    rxnList1 = model.rxns;
    differentSetsFlag = false;
else
    if (isempty(rxnList1))
        rxnList1 = model.rxns;
    end
end

if (nargin < 4)
    rxnList2 = rxnList1;
    differentSetsFlag = false;
else
    if (isempty(rxnList2))
        rxnList2 = rxnList1;
        differentSetsFlag = false;
    else
        differentSetsFlag = true;
    end
end
if (nargin < 5)
    verbFlag = false;
end

nrxns = length(model.rxns);

% Run single rxn deletions first to figure out lethal rxns
fprintf('Single deletion analysis to remove lethal rxns\n');
[singleRatio1,singleRate1,grRateWT] = singleRxnDeletion(model,method,rxnList1,verbFlag);
singleLethal1 = (singleRatio1 < 1e-9);
rxnListOrig1 = rxnList1;
rxnListOrig2 = rxnList1;
rxnList1 = rxnList1(~singleLethal1);
singleRate = singleRate1(~singleLethal1);
[tmp,listMap1] = ismember(rxnListOrig1,rxnList1);
fprintf('%d non-lethal rxns\n',length(rxnList1));

% Repeat the analysis for the second set of rxns
if (differentSetsFlag)
    fprintf('Single deletion analysis to remove lethal rxns from rxn set 2\n');
    [singleRatio2,singleRate2,grRateWT] = singleRxnDeletion(model,method,rxnList2,verbFlag);
    singleLethal2 = (singleRatio2 < 1e-9);
    rxnListOrig2 = rxnList2;
    rxnList2 = rxnList2(~singleLethal2);
    [tmp,listMap2] = ismember(rxnListOrig2,rxnList2);
    fprintf('%d non-lethal rxns\n',length(rxnList2));
else
    rxnList2 = rxnList1;
    listMap2 = listMap1;
end

nDelrxns1 = length(rxnList1);
nDelrxns2 = length(rxnList2);

grRateKO = ones(nDelrxns1,nDelrxns2)*grRateWT;
grRatioDble = ones(nDelrxns1,nDelrxns2);

if (differentSetsFlag)
    nTotalPairs = nDelrxns1*nDelrxns2;
else
    nTotalPairs = nDelrxns1*(nDelrxns1-1)/2;
end

% Run double deletion analysis
delCounter = 0;
fprintf('Double rxn deletion analysis\n');
fprintf('Total of %d pairs to analyze\n',nTotalPairs);
showprogress(0,'Double rxn deletion analysis in progress ...');
t = cputime;
fprintf('Perc complete\tCPU time\n');
for rxnNo1 = 1:nDelrxns1
   
    % Find rxn index
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
            save doublerxnDeletionTmp.mat grRateKO
        end
        [~ ,rxnID2] = ismember(rxnList2{rxnNo2},model.rxns);
        % Use FBA/MOMA/lMOMA to calculate deletion strain growth rate
        constrRxnInd = [rxnID1; rxnID2];
        modelTmp = model;
        modelTmp.lb(constrRxnInd) = 0;
        modelTmp.ub(constrRxnInd) = 0;
        % Get double deletion growth rate
        switch method
            case 'lMOMA'
                solKO = linearMOMA(model,modelTmp,'max');
            case 'MOMA'
                solKO = MOMA(model,modelTmp,'max',false,true);
            otherwise
                solKO = solveFBAmodelCplex(modelTmp);
        end
        %solKO = optimizeCbModel(modelTmp,'max');
        if (solKO.stat > 0)
            grRateKO(rxnNo1,rxnNo2) = solKO.f;
            grRateKO(rxnNo2,rxnNo1) = solKO.f;
        else
            grRateKO(rxnNo1,rxnNo2) = 0;
            grRateKO(rxnNo2,rxnNo1) = 0;
        end
        if (verbFlag)
            fprintf('%4d\t%4.1f\t%10s\t%10s\t%9.3f\t%9.3f\n',delCounter,100*delCounter/nTotalPairs,rxnList1{rxnNo1},...
                rxnList2{rxnNo2},grRateKO(rxnNo1,rxnNo2),grRateKO(rxnNo1,rxnNo2)/grRateWT*100);
        end
        if (differentSetsFlag)
            grRateKO(rxnNo2,rxnNo1) = grRateKO(rxnNo1,rxnNo2);
        end
    end
end

% Reconstruct the entire matrix
for i = 1:length(rxnListOrig1)
    for j = 1:length(rxnListOrig2)
        if (listMap1(i) > 0 & listMap2(j) > 0)
            allGrRateKO(i,j) = grRateKO(listMap1(i),listMap2(j));
        else
            allGrRateKO(i,j) = 0;
        end
    end
end

grRatioDble = allGrRateKO/grRateWT;

grRateKO = allGrRateKO;
