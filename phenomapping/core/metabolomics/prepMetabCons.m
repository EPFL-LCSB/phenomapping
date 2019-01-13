function [NewModel, LCcons, tagFeas] = prepMetabCons(model, metNames, ...
    LCmin, LCmax)
% Integrates metabolomics data into a model (with TFA structure as defined 
% in matTFA): one provides a metabolite name and its concentration value
% (mol/Lcell). The constraint for the metabolite concentration will be 
% integrated for that metabolite in all intracellular compartments where
% it appears.
%
% USAGE:
%
%    [NewModel, LCcons, tagFeas] = prepMetabCons(model, metNames, LCmin, LCmax)
%
% INPUT:
%    model:           TFA model structure
%    metNames:        List of metabolite names for which metabolomics data
%                     is available. The names should match the metNames in
%                     the model
%    LCmin:           LN(lowest concentrations values: normally, mean-SD)
%                     for the list of metabolites in metNames
%    LCmax:           LN(highest concentrations values: normally, mean+SD)
%                     for the list of metabolites in metNames
%
% OUTPUTS:
%    NewModel:        Model with metabolomics data integrated: for the
%                     metNames that were found in the model and in all
%                     compartments where this metabolite appears.
%    LCcons:          Summary of the data integrated in the standard format
%                     used to integrate metabolomics data into a TFA model
%                     LCcons = [LC_metID, LN(min concentration), LN(max
%                     concentration)]
%    tagFeas:         Tag to define whether the model is feasible upon 
%                     integration of metabolomics data.
%
% .. Author:
% Anush Chiappino-Pepe 2015
%

if ~any(ismember(metNames,model.metNames))
    error('none of the metabolites in the metabolomic data set is part of the model!')
elseif all(ismember(metNames,model.metNames))
    fprintf('all metabolites in the metabolomic data set are part of the model\n');
else
    warning('not all metabolites in the metabolomic data set are part of the model\n');
    metNames = metNames(ismember(metNames,model.metNames));
    LCmin = LCmin(ismember(metNames,model.metNames));
    LCmax = LCmax(ismember(metNames,model.metNames));
end

% rounding to 4th decimal
LCmin = roundn(LCmin,-4);
LCmax = roundn(LCmax,-4);

if ~isfield(model,'metCompSymbol')
    warning('the model does not contain the field metCompSymbol. One will be created. It is recommended to run initTestPhenoMappingModel first!')
    model = addMetCompSymbol(model);
end

% find metabolites in all intracellular compartments and assign
% concentration values
constraints = [];
LC = [];
for i = 1:length(metNames)
    tmp = find(ismember(model.metNames,metNames(i))); % find metabolite name
    tmp = tmp(~ismember(model.metCompSymbol(tmp),'e')); % discard extracellular metabolite
    constraints = [constraints; model.mets(tmp)]; % find the metabolite ID
    
    LCi = [];
    LCi(1:length(model.mets(tmp)),1) = LCmin(i); % assign to the metabolite the concentration constraint in all the intracellular compartments
    LCi(1:length(model.mets(tmp)),2) = LCmax(i);
    LC = [LC; LCi];
end

% defining standard format of metabolomic sets for matTFA
LCcons = [strcat('LC_',constraints), num2cell(LC)];
% integrating metabolomics set into the model
NewModel = loadConstraints(model, LCcons);
sol = optimizeThermoModel(NewModel);

tagFeas = 1;
if isnan(sol.val) || isempty(sol.val) || (sol.val<1E-3)
    tagFeas = 0;
    warning('You need to generate a reduced metabolomics data set for integration. Check tutorial_issues.m');
end