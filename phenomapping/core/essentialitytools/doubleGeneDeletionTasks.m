function [grRatio, grRateKO, grRateWT, delRxns, impactTasks] = doubleGeneDeletionTasks(model, method, gpairsList, verbFlag, flagTasks, essThr)
% Performs single gene deletion analysis using FBA, MOMA or linearMOMA
%
% USAGE:
%
%    [grRatio, grRateKO, grRateWT, delRxns, impactTasks] = doubleGeneDeletionTasks(model, method, gpairsList, verbFlag, flagTasks, essThr)
%
% INPUT:
%    model:           COBRA model structure including gene-reaction associations
%
% OPTIONAL INPUTS:
%    method:          Either 'FBA', 'MOMA' or 'lMOMA' (Default = 'FBA')
%    gpairsList:      Pairs of genes to be deleted (default = none)
%    verbFlag:        Verbose output (Default false)
%    flagTasks:       Determine which BBBs cannot be produced upon knockout
%                     (default = false)
%    essThr:          Threshold on growth below which the KO is considered
%                     to be lethal; required for flagTasks (default = 0.1)
%
%
% OUTPUTS:
%    grRatio:         Computed growth rate ratio between deletion strain and wild type
%    grRateKO:        Deletion strain growth rates (1/h)
%    grRateWT:        Wild type growth rate (1/h)
%    delRxns:         List of deleted reactions for each gene `KO`
%    impactTasks:     list of BBBs and its production upon KO of gpairsList
%
% .. Author:
%       - Markus Herrgard 8/7/06
%       - Aurich/Thiele 11/2015 unique gene deletion option (delete all alternate transcripts and if solKO.stat not 1 or 5, grRateKO(i) = NaN;)
%       - Karthik Raman 06/2017 speeding up gene deletion based on github.com/RamanLab/FastSL
%       - Anush Chiappino-Pepe 08/2017 check tasks and doubleKO

if (nargin < 2)
    method = 'FBA';
end
if (nargin < 3)
    gpairsList = {};
end
if (nargin < 4)
    verbFlag = false;
end
if (nargin < 5)
    flagTasks = false;
end
if (nargin < 6)
    essThr = 0.1;
end

if (isempty(gpairsList))
    error('provide a list of SL to test tasks');
end

nDelGenes = size(gpairsList,1);
impactTasks=cell(nDelGenes,2);
impactTasks(:,1) = strcat(gpairsList(:,1),'|',gpairsList(:,2));

if strcmp(method,'TFA')
    solWT = optimizeThermoModel(model);
    grRateWT = solWT.val;
else
    solWT = solveFBAmodelCplex(model);
    grRateWT = solWT.f;
end

grRateKO = ones(nDelGenes,1)*grRateWT;
% Assign the WT flux distribution to all deletions; those that differ
% will be replaced in the loop below
delRxns = cell(nDelGenes,1);
if (verbFlag)
    fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Growth rate','Rel. GR');
end
showprogress(0,'Single gene deletion analysis in progress ...');
for i = 1:nDelGenes
    showprogress(i/nDelGenes);
    if strcmp(method,'TFA')
        [modelDel,~,constrRxnNames1] = thermoDeleteModelGenes(model,gpairsList{i,1});
        [modelDel,~,constrRxnNames2] = thermoDeleteModelGenes(modelDel,gpairsList{i,2});
    else
        [modelDel,~,constrRxnNames1] = deleteModelGenes(model,gpairsList{i,1});
        [modelDel,~,constrRxnNames2] = deleteModelGenes(modelDel,gpairsList{i,2});
    end
    delRxns{i} = [constrRxnNames1;constrRxnNames2];
    % If all the reactions being deleted carried no flux in WT,
    % deleting them cannot affect the flux solution.
    switch method
        case 'lMOMA'
            solKO = linearMOMA(model,modelDel,'max');
        case 'MOMA'
            solKO = MOMA(model,modelDel,'max',false,true);
        case 'TFA'
            solKO = optimizeThermoModel(modelDel);
        otherwise
            %                     solKO = optimizeCbModel(modelDel,'max');
            solKO = solveFBAmodelCplex(modelDel);
    end
    if (solKO.stat == 1)
        if strcmp(method,'TFA')
            grRateKO(i) = solKO.val;
        else
            grRateKO(i) = solKO.f;
        end
        if flagTasks
            if grRateKO(i) < essThr*grRateWT
                [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,grRateWT);
            else
                impactTasks{i,2} = [];
            end
        end
    else
        grRateKO(i) = NaN;
        if flagTasks
            [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
        end
    end
    if (verbFlag)
        fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/nDelGenes,gpairsList{i},grRateKO(i),grRateKO(i)/grRateWT*100);
    end
end

grRatio = grRateKO/grRateWT;
end
