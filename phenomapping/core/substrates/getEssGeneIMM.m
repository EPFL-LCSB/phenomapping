function [essIMM, solOpt, models] = getEssGeneIMM(model, DPs, ...
    modelmilp, method, essThr, essGenesFBA, flagUpt, jointIMM, filename)
% Identifies essential genes in IMM/IMS
%
% USAGE:
%
%    [essIMM, solOpt, models] = getEssGeneIMM(model, DPs, modelmilp, method, essThr, essGenesFBA, flagUpt, jointIMM, filename)
%
% INPUT:
%    model:           TFA model structure (the "model" input to 
%                     analysisIMM.m)
%    DPs:             Directionality profile matrix with alternatives in
%                     each column
%    modelmilp:       TFA model with MILP structure for IMM/IMS analysis
%
% OPTIONAL INPUTS:
%    method:          Either 'FBA' or 'TFA' (default = 'FBA')
%    essThr:          Threshold on growth below which the KO is considered
%                     to be lethal; required for flagTasks (default = 0.1)
%    essGenesFBA:     For TFA gene essentiality, the list of essential
%                     genes from FBA for each IMM - output of this function
%                     (default = empty)
%    flagUpt:         True to identify the In silico Minimal Media (IMM).
%                     Else it gets the In silico Minimal Secretion (IMS).
%                     (default = true)
%    jointIMM:        Medium composition of jointIMM or any other at which 
%                     you might want to test gene essentiality (default =
%                     empty - to test essentiality at the IMMs of the DPs)
%    filename:        Name used to save DPs (default = 'imm')
%
% OUTPUTS:
%    essIMM:          Essentiality per IMM
%    solOpt:          Growth achieved in each IMM (all should be > 0)
%    models:          Models with constraints from each IMM integrated
%
% .. Author:
% Anush Chiappino-Pepe 2015
%

if (nargin < 4)
    method = 'FBA';
end
if (nargin < 5)
    essThr = 0.1;
end
if (nargin < 6)
    essGenesFBA = [];
end
if (nargin < 7)
    flagUpt = 1; %analysis of IMM
end
if (nargin < 8)
    jointIMM = [];
end
if (nargin < 9)
    filename = 'PhenoMappingSubstrates';
end

gr = optimizeThermoModel(model);
gr = gr.val;

% Extract the information from the DPs to know which metabolites were uptaken (use=0) and which secreted (use=1)
if ~isempty(jointIMM)
    [~, Mets] = extractInfoIMMDPs(modelmilp, DPs, modelmilp.indUSE);
    MatrixInfo = ones(size(DPs(modelmilp.indUSE,1),1),1);
    MatrixInfo(ismember(Mets(:,1),jointIMM)) = 0;
else
    MatrixInfo = DPs(modelmilp.indUSE,:);
end
Metsdrains = strrep(modelmilp.varNames(modelmilp.indUSE),'BFUSE_R_','');
Metsdrains = strrep(Metsdrains,'BFUSE_F_','');

models = cell(size(MatrixInfo,2),1);
solOpt = cell(size(MatrixInfo,2),1);
essIMM = cell(size(MatrixInfo,2),1);

% Essentiality analysis per DP
for i = 1:size(MatrixInfo,2)
    models{i,1} = model;
    if flagUpt
        % Define lb=0 if use=1 (secretion) for the exchange reactions if it is IMM
        models{i,1}.lb(ismember(models{i,1}.rxns,Metsdrains(MatrixInfo(:,i)>0.9))) = 0;
        models{i,1}.var_lb(ismember(models{i,1}.varNames,strcat('NF_',Metsdrains(MatrixInfo(:,i)>0.9)))) = 0;
    else
        models{i,1}.ub(ismember(models{i,1}.rxns,Metsdrains(MatrixInfo(:,i)>0.9))) = 0;
        models{i,1}.var_ub(ismember(models{i,1}.varNames,strcat('NF_',Metsdrains(MatrixInfo(:,i)>0.9)))) = 0;
    end
    % Essentiality studies for each model
    if strcmp(method,'FBA')
        changeCobraSolver('cplex_direct','LP');
        [~, grRateKO] = singleGeneDeletionTasks(models{i,1}, method, models{i,1}.genes);
        grRateKO(isnan(grRateKO)) = 0;
        solOpt{i,1} = solveFBAmodelCplex(models{i,1});
        solOpt{i,1} = solOpt{i,1}.f;
        essIMM{i,1} = models{i,1}.genes(grRateKO < essThr*gr);
    elseif strcmp(method,'TFA')
        geneList = models{i,1}.genes(~ismember(models{i,1}.genes,essGenesFBA{i,1}));
        [~, grRateKO] = thermoSingleGeneDeletion(models{i,1}, method, geneList, 0, 0, 0, 0.1, model.indNF);
        grRateKO(isnan(grRateKO)) = 0;
        solOpt{i,1} = optimizeThermoModel(models{i,1});
        solOpt{i,1} = solOpt{i,1}.val;
        essIMM{i,1} = geneList(grRateKO < essThr*gr);
    else
        error('only FBA and TFA have been implemented')
    end
    if rem(i,50) == 0
        save(strcat(filename,'_ess.mat'), 'essIMM');
    end
end
