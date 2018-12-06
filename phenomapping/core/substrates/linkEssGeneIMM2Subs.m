function [subsToGenes, essIMMaddToIRM] = ...
    linkEssGeneIMM2Subs(modelmilp, essIMM, DPs, model, essThr, essIRM, ...
    NumAlt, time, jointIMM, filename)
% Links essentiality of a gene with the subsToGenes missing in the IMM
%
% USAGE:
%
%    [subsToGenes, essIMMaddToIRM, essSubsToAdd] = linkEssGeneIMM2Subs(modelmilp, essIMM, DPs, model, essThr, essIRM, NumAlt, time, jointIMM, filename)
%
%
% INPUT:
%    modelmilp:       TFA model with MILP structure for IMM/IMS analysis
%    essIMM:          Essentiality per IMM as output of getEssentialityIMM.m                  
%    DPs:             Directionality profile matrix with alternatives in
%                     each column
%    model:           TFA model structure (the "model" input to 
%                     analysisIMM.m)
%
%
% OPTIONAL INPUTS:
%    essThr:          Threshold on growth below which the KO is considered
%                     to be lethal; required for essIRM (default = 0.1)
%    essIRM:          Essentiality in rich medium (default = tested with
%                     thermoSingleGeneDeletion.m)
%    NumAlt:          Number of alternatives to generate (default = 1)
%    time:            Time in seconds (default = empty / not time limit)
%    jointIMM:        Medium composition of jointIMM or any other at which 
%                     you might want to test gene essentiality (default =
%                     empty - to test essentiality at the IMMs of the DPs)
%    filename:        Name used to save DPs (default = 'imm')
%
% OUTPUTS:
%    subsToGenes:     Substrates whose abscence in the medium determine
%                     the essential function of the genes in essIMMaddToIRM
%    essIMMaddToIRM:  Essential genes in IMM on the top of essential genes
%                     in rich medium
%
% NOTE: model, essIRM and essIMM should be given for the same model and
% analysis method (FBA, TFA or TFA with metabolomics)
%
% .. Author:
% Anush Chiappino-Pepe 2015
%

gr = optimizeThermoModel(model);
gr = gr.val;

if (nargin < 5)
    essThr = 0.1;
end
if (nargin < 6)
    fprintf('TFA essentiality will be performed\n');
    grRatioGeneTFA = thermoSingleGeneDeletion(model, 'TFA', model.genes, 0, 0, 0, essThr, model.indNF);
    grRatioGeneTFA(isnan(grRatioGeneTFA)) = 0; %by default NaN is considered an essential KO
    essIRM = model.genes(grRatioGeneTFA<essThr);
end
if (nargin < 7)
    NumAlt = 10;
end
if (nargin < 8)
    time = [];
end
if (nargin < 9)
    jointIMM = [];
end
if (nargin < 10)
    filename = 'PhenoMappingSubstrates';
end

% pre-allocate outputs
subsToGenes = cell(size(essIMM,1),1);
essIMMaddToIRM = cell(size(essIMM,1),1);

% get indexes of integer and continuous variables of the milp
indUSE = modelmilp.indUSE;
[~, indCont] = ismember(strrep(modelmilp.varNames(indUSE),'BFUSE_',''),modelmilp.varNames); % get indexes of continuous variable associated to integers
% extract part of matrix with milp solution from IMM
if ~isempty(jointIMM)
    % this is for the analysis of the jointIMM or any other media of
    % interest
    [~, Mets] = extractInfoIMMDPs(modelmilp, DPs, indUSE);
    MatrixInfo = ones(size(DPs(indUSE,1),1),1);
    MatrixInfo(ismember(Mets(:,1),jointIMM)) = 0;
else
    MatrixInfo = DPs(indUSE,:);
end

% set lb in biomass if this was not defined before
if modelmilp.var_lb(ismember(modelmilp.varNames,strcat('F_',modelmilp.rxns(modelmilp.c==1))))<1e-7
    warning('there is no defined lb in the objective function')
    modelmilp.var_lb(ismember(modelmilp.varNames,strcat('F_',modelmilp.rxns(modelmilp.c==1))))=1e-3;
end

% check there is no cut constraint limiting the length of the milp solution
if ismember('CUT_0',modelmilp.constraintNames)
    warning('the model you have sent to this function has CUT constraints integrated')
    modelmilp.rhs(ismember(modelmilp.constraintNames,'CUT_0'))=0;
end

% get pool of all genes essential in the IMMs in addition to IRM
[~, commonIMM, notCommonIMM] = getOverlapSets(essIMM);
jointAddIMM = unique([commonIMM; notCommonIMM]);
jointAddIMM = jointAddIMM(~ismember(jointAddIMM,essIRM));
% we will get track of the genes analyzed to avoid repeating analysis
taggene = zeros(length(jointAddIMM),1);

for z = 1:size(MatrixInfo,2)
    % find the additional essential genes of each IMM wrt the IRM 
    ess_add = essIMM{z,1}(~ismember(essIMM{z,1},essIRM));
    
    % avoid analyzing the same gene twice: taggene=1 was done already
    if (z >= 2) && ~isempty(jointAddIMM)
        ess_add = ess_add(ismember(ess_add, jointAddIMM(taggene==0)));
    end
    essIMMaddToIRM{z,1} = ess_add;
    save(strcat(filename,'_addGenes.mat'), 'essIMMaddToIRM');
    
    % get the subset of integer/continuous variables to study (the
    % substrates that are not part of the IMM/IMS)
    intUSE = indUSE(MatrixInfo(:,z)>0.9);
    intCont = indCont(MatrixInfo(:,z)>0.9);
    
    for j = 1:length(ess_add)
        if ~isempty(jointAddIMM)
            % tag with 1 the genes analyzed
            taggene(ismember(jointAddIMM,ess_add(j))) = 1;
        end
        % block essential gene in this IMM
        modelKO = thermoDeleteModelGenes(modelmilp,ess_add(j));
        % initialize the integer variables for the milp
        modelKO.f = zeros(length(modelKO.varNames),1);
        
        % verify that gene is essential at this IMM
        tt = modelKO;
        tt.f(ismember(tt.varNames,strcat('F_',tt.rxns(tt.c==1)))) = 1;
        tt.var_ub(intCont) = 0; % block uptakes or secretions active in this IMM
        ss = optimizeThermoModel(tt);
        
        if (isempty(ss.val) || ss.val<essThr*gr || isnan(ss.val)) % if essential continue to identify substrates that allow survival without this gene
            % define integers only for the substrates that are not part of
            % the IMM of study!
            modelKO.f(intUSE) = 1;
            modelKO.objtype = -1; %minimize:1, maximize:-1 (just for verification!)
            cDPs = findDPMax(modelKO, NumAlt, intUSE, time, 0, filename); % optimize!
            
            if (isempty(cDPs))
                subsToGenes{z,1}{j,1} = 'problem: should be essential in the wtmodel';
            end
            for i = 1:size(cDPs,2)
                rxns = modelmilp.varNames(intUSE(cDPs(intUSE,i)<0.1));
                if (~isempty(rxns))
                    rxnsEss=strrep(rxns, 'BFUSE_R_', '');
                    rxnsEss=strrep(rxnsEss, 'BFUSE_F_', '');
                    subsToGenes{z,1}{j,i} = printRxnFormula(modelmilp,rxnsEss,0,0,1);
                else
                    subsToGenes{z,1}{j,i} = '';
                end
            end
            save(strcat(filename,'_subsToGenes.mat'), 'subsToGenes');
        else
            warning('check out! a gene marked as essential is not essential in the IMM')
        end
    end
end