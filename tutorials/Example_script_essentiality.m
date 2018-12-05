%% Example script for essentiality analysis of genes, reactions, and metabolites with FBA, MOMA, TFA and tMOMA


%% Inputs
% obligatory inputs
tmodel = human;
essThr = 0.1; % essentiality threshold
flagTasks = 0; % true if you want to perform essentiality per task

% optional inputs
description = 'human';

% select solver for FBA and MOMA (default cplex)
changeCobraSolver('cplex_direct','LP');
changeCobraSolver('ibm_cplex','QP');


%% 
if ~isfield(tmodel,'rules')
    tmodel = generateRules(tmodel);
end

sol = solveFBAmodelCplex(tmodel);
solt = optimizeThermoModel(tmodel);

indNF = getAllVar(tmodel,{'NF'});
if isempty(indNF)
    tmodel = addNetFluxVariables(tmodel);
    indNF = getAllVar(tmodel,{'NF'});
end

%% RxnKO with FBA
[grRatio_rxn,grRateKO_rxn] = singleRxnDeletionTasks(tmodel,'FBA',tmodel.rxns);

if any(isnan(grRateKO_rxn)) %keep track of the NaN KO by comparing essential_rxns_fbaNaN with essential_rxns_fba
    grRateKO_rxnNaN = grRateKO_rxn;
    essential_rxns_fbaNaN = tmodel.rxns(grRateKO_rxnNaN<essThr*sol.f);
    yesrxnfbaNaN = ismember(tmodel.rxns,essential_rxns_fbaNaN);
end

grRateKO_rxn(isnan(grRateKO_rxn)) = 0; %by default NaN is considered an essential KO
essential_rxns_fba = tmodel.rxns(grRateKO_rxn<essThr*sol.f);
yesrxnfba = ismember(tmodel.rxns,essential_rxns_fba);

%% GeneKO with FBA
[grRatio_gene,grRateKO_gene] = singleGeneDeletionTasks(tmodel,'FBA',tmodel.genes);

if any(isnan(grRateKO_gene)) %keep track of the NaN KO by comparing essential_genes_fbaNaN with essential_genes_fba
    grRateKO_geneNaN = grRateKO_gene;
    essential_genes_fbaNaN = tmodel.genes(grRateKO_geneNaN<essThr*sol.f);
    yesgenefbaNaN = ismember(tmodel.genes,essential_genes_fbaNaN);
end

grRateKO_gene(isnan(grRateKO_gene)) = 0; %by default NaN is considered an essential KO
essential_genes_fba = tmodel.genes(grRateKO_gene<essThr*sol.f);
yesgenefba = ismember(tmodel.genes,essential_genes_fba);

%% check tasks for essential rxns and genes with FBA
if flagTasks
    [~, ~, ~, ~, ~, ~, impactTasks_rxn] = singleRxnDeletionTasks(tmodel, 'FBA', essential_rxns_fba, 0, flagTasks, essThr);
    essTasks_rxn = extractTasksImpact(impactTasks_rxn(:,2), essential_rxns_fba);
    [~, ~, ~, ~, ~, ~, impactTasks_gene] = singleGeneDeletionTasks(tmodel, 'FBA', essential_genes_fba, 0, 0, flagTasks, essThr);
    essTasks_gene = extractTasksImpact(impactTasks_gene(:,2), essential_genes_fba);
end
%% RxnKO with TFA
[grRatio_rxntfa,grRateKO_rxntfa] = thermoSingleRxnDeletion(tmodel, 'TFA', tmodel.rxns, 0, 0, essThr, indNF);

if any(isnan(grRateKO_rxntfa)) %keep track of the NaN KO by comparing essential_rxns_tfaNaN with essential_rxns_tfa
    grRateKO_rxntfaNaN = grRateKO_rxntfa;
    essential_rxns_tfaNaN = tmodel.rxns(grRateKO_rxntfaNaN(:,1)<essThr*solt.val);
    yesrxntfaNaN = ismember(tmodel.rxns,essential_rxns_tfaNaN);
end

grRateKO_rxntfa(isnan(grRateKO_rxntfa)) = 0; %by default NaN is considered an essential KO
essential_rxns_tfa = tmodel.rxns(grRateKO_rxntfa(:,1)<essThr*solt.val); 
yesrxntfa = ismember(tmodel.rxns,essential_rxns_tfa);

%% GeneKO with TFA
[grRatio_genetfa,grRateKO_genetfa] = thermoSingleGeneDeletion(tmodel, 'TFA', tmodel.genes, 0, 0, 0, essThr, indNF);

if any(isnan(grRateKO_rxntfa)) %keep track of the NaN KO by comparing essential_genes_tfaNaN with essential_genes_tfa
    grRateKO_genetfaNaN = grRateKO_genetfa;
    essential_genes_tfaNaN = tmodel.genes(grRateKO_genetfaNaN(:,1)<essThr*solt.val);
    yesgenetfaNaN = ismember(tmodel.genes,essential_genes_tfaNaN);
end

grRateKO_genetfa(isnan(grRateKO_genetfa)) = 0; %by default NaN is considered an essential KO
essential_genes_tfa = tmodel.genes(grRateKO_genetfa(:,1)<essThr*solt.val);
yesgenetfa = ismember(tmodel.genes,essential_genes_tfa);

%% check tasks for essential rxns and genes with TFA
if flagTasks
    [~, ~, ~, ~, ~, ~, impactTasks_rxn_tfa] = thermoSingleRxnDeletion(tmodel, 'TFA', essential_rxns_tfa, 0, flagTasks, essThr, indNF);
    essTasks_rxn_tfa = extractTasksImpact(impactTasks_rxn_tfa(:,2), essential_rxns_tfa);
    [~, ~, ~, ~, ~, ~, impactTasks_gene_tfa] = thermoSingleGeneDeletion(tmodel, 'TFA', essential_genes_tfa, 0, 0, flagTasks, essThr, indNF);
    essTasks_gene_tfa = extractTasksImpact(impactTasks_gene_tfa(:,2), essential_genes_tfa);
end

%% Double GeneKO with FBA
[grRatioDble,grRateKO] = doubleGeneDeletion(tmodel,'FBA',tmodel.genes(not(ismember(tmodel.genes,essential_genes_fba))));
genesDouble = tmodel.genes(not(ismember(tmodel.genes,essential_genes_fba)));
[pairsgenes, grRateEffectpairs] = extractDoubleGeneMatrix(tmodel,grRateKO, genesDouble);
EssentialPairs = pairsgenes(grRateEffectpairs(:,1)<essThr*sol.f,:);

%% Double GeneKO with TFA
% suggestion: send here a model (reducedModel) that does not contain blocked reactions and
% their associated genes
geneList_doublethermo = tmodel.genes(not(ismember(tmodel.genes,essential_genes_tfa)));
[tgrRateKOmatrix]=thermoDoubleGeneDeletion(tmodel,geneList_doublethermo,essThr);
[pairsgenes_thermo, grRateEffectpairs_thermo] = extractDoubleGeneMatrix(tmodel,tgrRateKOmatrix,geneList_doublethermo);

if any(isnan(grRateEffectpairs_thermo)) %keep track of the NaN KO by comparing EssentialPairs_thermoNaN with EssentialPairs_thermo
    grRateEffectpairs_thermoNaN = grRateEffectpairs_thermo;
    EssentialPairs_thermoNaN = pairsgenes_thermo(grRateEffectpairs_thermoNaN(:,1)<essThr*solt.val,:);
end

grRateEffectpairs_thermo(isnan(grRateEffectpairs_thermo)) = 0; %by default NaN is considered an essential KO
EssentialPairs_thermo = pairsgenes_thermo(grRateEffectpairs_thermo(:,1)<essThr*solt.val,:);

%% RxnKO with MOMA
[grRatio_rxnmm,grRateKO_rxnmm] = singleRxnDeletion(tmodel,'MOMA');
essential_rxns_moma = tmodel.rxns(grRateKO_rxnmm<essThr*sol.f);
yesrxnMOMA = ismember(tmodel.rxns,essential_rxns_moma);

%% GeneKO with MOMA
[grRatio_genemm,grRateKO_genemm] = singleGeneDeletion(tmodel,'MOMA');
essential_genes_moma = tmodel.genes(grRateKO_genemm<essThr*sol.f);
yesgeneMOMA = ismember(tmodel.genes,essential_genes_moma);

%% GeneKO with tMOMA
[grRatio_genetfa_moma,grRateKO_genetfa_moma] = thermoSingleGeneDeletion(tmodel, 'tMOMA');
essential_genes_tmoma = tmodel.genes(grRateKO_genetfa_moma(:,1)<essThr*solt.val);
yesgenetMOMA = ismember(tmodel.genes,essential_genes_tmoma);

%% MetKO with FBA
[grRateKO_metfba, grRatioKO_metfba] = testMetEss(tmodel, tmodel.mets);

if any(isnan(grRateKO_metfba)) %keep track of the NaN KO by comparing essential_met_fbaNaN with essential_mets_fba
    grRateKO_metfbaNaN = grRateKO_metfba;
    essential_met_fbaNaN = tmodel.mets(grRateKO_metfbaNaN<essThr*sol.f);
    yesmetfbaNaN = ismember(tmodel.mets,essential_mets_fbaNaN);
end

grRateKO_metfba(isnan(grRateKO_metfba)) = 0; %by default NaN is considered an essential KO
essential_mets_fba = tmodel.mets(grRateKO_metfba<essThr*sol.f);
yesmetfba = ismember(tmodel.mets,essential_mets_fba);

%% MetKO with TFA
[grRateKO_mettfa, grRatioKO_mettfa] = testMetEssThermo(tmodel, tmodel.mets);

if any(isnan(grRateKO_mettfa)) %keep track of the NaN KO by comparing essential_met_tfaNaN with essential_mets_tfa
    grRateKO_mettfaNaN = grRateKO_mettfa;
    essential_met_tfaNaN = tmodel.mets(grRateKO_mettfaNaN<essThr*sol.f);
    yesmettfaNaN = ismember(tmodel.mets,essential_met_tfaNaN);
end

grRateKO_mettfa(isnan(grRateKO_mettfa)) = 0; %by default NaN is considered an essential KO
essential_mets_tfa = tmodel.mets(grRateKO_mettfa<essThr*sol.f);
yesmettfa = ismember(tmodel.mets,essential_mets_tfa);

