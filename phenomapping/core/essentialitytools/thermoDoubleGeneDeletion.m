function [grRateKOmatrix, geneList] = thermoDoubleGeneDeletion(model, geneList, essThr)
% Performs double gene deletion analysis using TFA
%
% USAGE:
%
%    [grRateKOmatrix] = thermoDoubleGeneDeletion(model, geneList, essThr)
%
% INPUT:
%    model:           TFA model structure including gene-reaction associations
%
% OPTIONAL INPUTS:
%    geneList:        List of genes to be deleted (default = all genes)
%    essThr:          Essentiality threshold for single gene deletion
%
% OUTPUTS:
%    grRateKOmatrix:  Matrix with the growth value upon double knockout of
%                     genes in geneList
%    geneList:        Genes tested for doubleKO
%
% .. Author:
%       - Anush Chiappino-Pepe 31/8/17
%
% Note: analyze the grRateKOmatrix with extractDoubleGeneMatrix.m

if nargin < 2
    geneList = [];
end
if nargin < 3
    essThr = 0.1;
end

sol = optimizeThermoModel(model);
indNF = getAllVar(model,{'NF'});
if isempty(indNF)
    model = addNetFluxVariablesNEW(model);
end

% if not provided, calculate essential genes with tfba (single KO)
% Extract the genes that are not essential with single thermo gene knockout
if isempty(geneList)
    [~, grRateKO_genetfa] = thermoSingleGeneDeletion(tmodel, 'TFA', model.genes);
    essential_genes_tfba = tmodel.genes(grRateKO_genetfa(:,1) < essThr*sol.val);
    geneList = model.genes(not(ismember(model.genes, essential_genes_tfba)));
end

delRxns = cell(length(geneList),1);
for i = 1:length(geneList)
    [~, ~, constrRxnNames] = thermoDeleteModelGenes(model, geneList{i});
    delRxns{i} = constrRxnNames;
end

% grRateKOmatrix will be a simetric matrix; avoid calculating twice the
% same by discarding elements = 1000 (lower triangular matrix)
grRateKOmatrix = tril(1000*ones(length(geneList),length(geneList)));

% For each of these genes create a model without it and perform single gene
% deletion
% Study the combined effect

showprogress(0,'Double gene knockout analysis in progress ...');
for i = 1:length(geneList)
    disp(i/length(geneList));
    [~ , geneInd1] = ismember(geneList(i),model.genes);
    % Identify Jzrxns, the set of reactions that do not carry a flux in TFA
    % minnorm with this gene knocked out; none of these can be lethal
    [modelDel] = thermoDeleteModelGenes(model,geneList{i});
%     solWTtfa = optimizeThermoModel(modelDel,true);
%     Jzrxns = model.rxns(solWTtfa.x(indNF)==0);

    
    for j = 1:length(geneList)
        showprogress(j/length(geneList));
        if grRateKOmatrix(i,j)==1000
        else
            [~ , geneInd2] = ismember(geneList(j),model.genes);
            rxnInd = find(any(model.rxnGeneMat(:,[geneInd1 geneInd2]),2));
            if (~isempty(rxnInd))
                % Initialize the state of all genes
                x = true(length(model.genes),1);
                x(geneInd1) = false;
                x(geneInd2) = false;
                constrainRxn = false(length(rxnInd),1);
                % Figure out if any of the reaction states is changed
                for rxnNo = 1:length(rxnInd)
                    if (~eval(model.rules{rxnInd(rxnNo)}))
                        constrainRxn(rxnNo) = true;
                    end
                end
                if (any(constrainRxn)) %&& ~all(ismember(delRxns{j},Jzrxns)))
                    constrRxnInd = rxnInd(constrainRxn);
                    [~ , ba] = ismember(strcat('NF_', model.rxns(constrRxnInd)), model.varNames);
                    dmodel = model;
                    dmodel.var_ub(ba) = 0;
                    dmodel.var_lb(ba) = 0;
                    soli = optimizeThermoModel(dmodel);
                    if sol.stat == 1
                        grRateKOmatrix(i,j) = soli.val;
                    else
                        grRateKOmatrix(i,j) = NaN;
                    end
                else
                    grRateKOmatrix(i,j) = sol.val;
                end
            else
                grRateKOmatrix(i,j) = sol.val;
            end
        end
    end
    save grRateKOmatrix grRateKOmatrix
end
