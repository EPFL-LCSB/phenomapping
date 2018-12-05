function [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, impactTasks] = singleGeneDeletionTasks(model, method, geneList, verbFlag, uniqueGene, flagTasks, essThr)
% Performs single gene deletion analysis using FBA, MOMA or linearMOMA
%
% USAGE:
%
%    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, impactTasks] = singleGeneDeletion(model, method, geneList, verbFlag, uniqueGene, flagTasks, essThr)
%
% INPUT:
%    model:           COBRA model structure including gene-reaction associations
%
% OPTIONAL INPUTS:
%    method:          Either 'FBA', 'MOMA' or 'lMOMA' (Default = 'FBA')
%    geneList:        List of genes to be deleted (default = all genes)
%    verbFlag:        Verbose output (Default false)
%    uniqueGene:      Run unique gene deletion (default = 0)
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
%    hasEffect:       Does a gene deletion affect anything (i.e. are any reactions
%                     removed from the model)
%    delRxns:         List of deleted reactions for each gene `KO`
%    fluxSolution:    FBA/MOMA/lMOMA fluxes for `KO` strains
%    impactTasks:     list of BBBs and its production upon KO of geneList
%
% .. Author:
%       - Markus Herrgard 8/7/06
%       - Aurich/Thiele 11/2015 unique gene deletion option (delete all alternate transcripts and if solKO.stat not 1 or 5, grRateKO(i) = NaN;)
%       - Karthik Raman 06/2017 speeding up gene deletion based on github.com/RamanLab/FastSL
%       - Anush Chiappino-Pepe 08/2017 check tasks
if (nargin < 2)
    method = 'FBA';
end
if (nargin < 3)
    geneList = model.genes;
else
    if (isempty(geneList))
        geneList = model.genes;
    end
end
if (nargin < 4)
    verbFlag = false;
end
if (nargin < 5)
    uniqueGene = 0;
end
if (nargin < 6)
    flagTasks = false;
end
if (nargin < 7)
    essThr = 0.1;
end

if (uniqueGene == 1)
    % detect whether there are alternate transcripts
    if ~isempty(strfind(model.genes{1},'.'))
        [geneList,rem] = strtok(model.genes,'.');
        geneList = unique(geneList);
        nGenes = length(geneList);
    else
        nGenes = length(model.genes);
    end
    nDelGenes = length(geneList);
    impactTasks=cell(nDelGenes,2);
    impactTasks(:,1) = geneList;
    
    solWT = solveFBAmodelCplex(model);
    grRateWT = solWT.f;
    
    grRateKO = ones(nDelGenes,1)*grRateWT;
    grRatio = ones(nDelGenes,1);
    hasEffect = true(nDelGenes,1);
    fluxSolution = zeros(length(model.rxns),nDelGenes);
    delRxns = cell(nDelGenes,1);
    if (verbFlag)
        fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Growth rate','Rel. GR');
    end
    showprogress(0,'Single gene deletion analysis in progress ...');
    for i = 1:nDelGenes
        showprogress(i/nDelGenes);
        if ~isempty(strfind(model.genes{1},'.'))
            % delete all alternate transcripts
            delGenes = model.genes(strmatch(geneList{i},model.genes));
            [modelDel,hasEffect(i),constrRxnNames] = deleteModelGenes(model,delGenes);
        else
            [modelDel,hasEffect(i),constrRxnNames] = deleteModelGenes(model,geneList{i});
        end
        delRxns{i} = constrRxnNames;
        if (hasEffect(i))
            switch method
                case 'lMOMA'
                    solKO = linearMOMA(model,modelDel,'max');
                case 'MOMA'
                    solKO = MOMA(model,modelDel,'max',false,true);
                otherwise
                    %                     solKO = optimizeCbModel(modelDel, 'max');
                    solKO = solveFBAmodelCplex(modelDel);
            end
            if (solKO.stat == 1 || solKO.stat == 5 )
                grRateKO(i) = solKO.f;
                fluxSolution(:,i) = solKO.x;
                if flagTasks
                    if grRateKO(i) < essThr*solWT.f
                        [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
                    else
                        impactTasks{i,2} = [];
                    end
                end
            else
                grRateKO(i) = NaN;
                fluxSolution(:,i) = nan(length(model.rxns),1);
                if flagTasks
                    [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
                end
            end
        end
        if (verbFlag)
            fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/nDelGenes,geneList{i},grRateKO(i),grRateKO(i)/grRateWT*100);
        end
    end


else
    nDelGenes = length(geneList);
    impactTasks=cell(nDelGenes,2);
    impactTasks(:,1) = geneList;
    
    solWT = solveFBAmodelCplex(model); %implement minNorm
    grRateWT = solWT.f;
    
    % Identify J_z, the set of reactions that do not carry a flux in solWT
    Jz = solWT.x==0;
    
    grRateKO = ones(nDelGenes,1)*grRateWT;
    hasEffect = true(nDelGenes,1);
    % Assign the WT flux distribution to all deletions; those that differ
    % will be replaced in the loop below
    fluxSolution = repmat(solWT.x, 1, nDelGenes);
    delRxns = cell(nDelGenes,1);
    if (verbFlag)
        fprintf('%4s\t%4s\t%10s\t%9s\t%9s\n','No','Perc','Name','Growth rate','Rel. GR');
    end
    showprogress(0,'Single gene deletion analysis in progress ...');
    for i = 1:nDelGenes
        showprogress(i/nDelGenes);
        [modelDel,hasEffect(i),constrRxnNames] = deleteModelGenes(model,geneList{i});
        delRxns{i} = constrRxnNames;
        % If all the reactions being deleted carried no flux in WT,
	    % deleting them cannot affect the flux solution. 
        if (hasEffect(i) && ~all(ismember(delRxns{i},model.rxns(Jz))))
            switch method
                case 'lMOMA'
                    solKO = linearMOMA(model,modelDel,'max');
                case 'MOMA'
                    solKO = MOMA(model,modelDel,'max',false,true);
                otherwise
                    %                     solKO = optimizeCbModel(modelDel,'max');
                    solKO = solveFBAmodelCplex(modelDel);
            end
            if (solKO.stat == 1)
                grRateKO(i) = solKO.f;
                fluxSolution(:,i) = solKO.x;
                if flagTasks
                    if grRateKO(i) < essThr*solWT.f
                        [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
                    else
                        impactTasks{i,2} = [];
                    end
                end
            else
                grRateKO(i) = NaN;
                fluxSolution(:,i) = nan(length(model.rxns),1);
                if flagTasks
                    [~, ~, impactTasks{i,2}] = checkBBBTasks(modelDel,method,essThr,[],0,0,solWT.f);
                end
            end
        end
        if (verbFlag)
            fprintf('%4d\t%4.0f\t%10s\t%9.3f\t%9.3f\n',i,100*i/nDelGenes,geneList{i},grRateKO(i),grRateKO(i)/grRateWT*100);
        end
    end
end

grRatio = grRateKO/grRateWT;
end
