function [model, checkList, tagReady] = initTestPhenoMappingModel(modeli, ...
    cplexPath, ReactionDBpath, tagThermo, rxnNoThermo, pathToSave)
% Initial check list to verify the model is ready for a PhenoMapping
% analysis
%
% USAGE:
%
%       [model, checkList, tagReady] = initTestPhenoMappingModel(modeli, cplexPath, ReactionDBpath, tagThermo, rxnNoThermo)
%
% INPUTS:
%    model:           model with FBA/TFA structure
%
% OPTIONAL INPUTS:
%    cplexPath:       Path to CPLEX (default = empty / to provide manually)
%    ReactionDBpath:  Path to database with all the thermodynamic 
%                     information from Alberty (default = empty / to
%                     provide manually)
%    tagThermo:       True if it is desired to convert the model to thermo
%                     with the matTFA toolbox (default = true)
%    rxnNoThermo:     List of rxns for which no thermo constraints should 
%                     be applied (default = empty / no reaction relaxed)
%
%
% OUTPUTS:
%    model:           model ready for PhenoMapping, unless it is missing
%                     major pieces of information and it requires a manual 
%                     curation
%    checkList:       list of tests done
%    tagReady:        True if model is ready for PhenoMapping
%
% Anush Chiappino-Pepe 2018
%

if (nargin < 2)
    cplexPath = [];
end
if (nargin < 3)
    ReactionDBpath = [];
end
if (nargin < 4)
    tagThermo = 1;
end
if (nargin < 5)
    rxnNoThermo = [];
end
if (nargin < 6)
    pathToSave = 'tmpresults/';
end

tagReady = 1;
model = modeli;

% PhenoMapping was developed to work with the solver CPLEX. We hence check 
% that you have CPLEX installed. In future releases, this repository
% will work with other solvers like gurobi.
fprintf('1: setting up cplex as the solver\n');
[solverOK,path_found] = addCplexPath(cplexPath);
if solverOK
    checkList{1} = strcat('ok: solver cplex path found ',path_found);
else
    checkList{1} = 'issue: solver cplex not found. Provide a valid path for cplex';
end

% test for feasibility of the model
sol = solveFBAmodelCplex(model);
if isnan(sol.f) || isempty(sol.f) || sol.f<1E-3
    error('the model provided is not feasible!')
end

% check the model structure required for thermo conversion
fprintf('2: checking model structure\n');
if ~isfield(model,'metCompSymbol')
    warning('the model does not contain the field metCompSymbol. The compartments will be searched automatically')
    model = addMetCompSymbol(model);
end

if tagThermo && (~isfield(model,'metSEEDID') || ~isfield(model,'CompartmentData'))
    warning('the model does not contain the fields metSEEDID or/and CompartmentData required in matTFA for conversion to thermo. Hence it wont be converted to thermo')
    tagThermo = 0;
end

% check that drains are all of the type A => (required for in silico
% minimal media analysis)
[model, flagChange] = putDrainsForward(model);
if flagChange
    checkList{2} = 'corrected: some drains in the model were not defined as A=> and the model was reconverted to a TFA structure';
else
    checkList{2} = 'ok: the model had all drains of the type A=>';
end

% check that the model has a TFA-friendly structure
if isfield(model,'A') && isfield(model,'f') && isfield(model,'var_lb') ...
        && isfield(model,'var_ub') && isfield(model,'rhs') ...
        && isfield(model,'constraintNames') && isfield(model,'varNames') ...
        && isfield(model,'constraintType') && isfield(model,'vartypes') ...
        && isfield(model,'objtype')
    checkList{3} = 'ok: the model already has a TFA structure';
    tagConv = 0;
else
    tagConv = 1;
    checkList{3} = 'corrected: the model was converted to a TFA structure';
end

% if the model does not have a TFA-friendly structure or it has drains like
% =>A create a TFA structure
if flagChange || tagConv
    if tagThermo
        % the TFA structure will have thermo constraints if requested and
        % possible (see first point: required fields)
        if isempty(ReactionDBpath)
            while isempty(ReactionDBpath)
                ReactionDBpath = input('Please provide the path to matTFA/thermoDatabases/thermo_data.mat and press enter\n... ','s');
            end
        end
        ReactionDB = load(ReactionDBpath);
        ReactionDB = ReactionDB.DB_AlbertyUpdate;
        modelt = prepModelforTFA(model, ReactionDB, model.CompartmentData);
        modelt = convToTFA(modelt, ReactionDB, rxnNoThermo);
        
        % it can happen that a model is not thermodynamically feasible because of
        % the directionality of a set of reactions. One should investigate these
        % cases in more detail. One option is to remove the thermodynamic
        % constraints for these reactions: this will be the input "rxnNoThermo".
        % Here, if the model is thermo infeasible, we will only apply FBA.
        sol = optimizeThermoModel(modelt);
        if isnan(sol.val) || isempty(sol.val) || sol.val<1E-3
            warning('the model is not thermodynamically feasible: we took out the thermo constraints. To work with thermo, please, investigate the reactions that make your model infeasible with the matTFA repository')
            model = convGEMToTFAStruc(model);
        end
    else
        model = convGEMToTFAStruc(model);
    end
end

sol = optimizeThermoModel(model);
if isnan(sol.val) || isempty(sol.val) || sol.val<1E-3
    error('the model is not feasible')
else
    fprintf('the model is feasible!\n');
end

% get the net fluxes associated with TFA structure
model.indNF = getAllVar(model,{'NF'});
if isempty(model.indNF)
    model = addNetFluxVariables(model);
    model.indNF = getAllVar(model,{'NF'});
end

% make sure the model has the structures required for gene essentiality
fprintf('3: checking fields for gene essentiality\n');
if isfield(model,'genes') && isfield(model,'grRules')
    checkList{4} = 'ok: the model has the fields for gene essentiality';
    if isfield(model,'rules')
        checkList{5} = 'ok: the field rules is present';
    else
        [model] = generateRules(model);
        checkList{5} = 'corrected: the field rules was added';
    end
else
    checkList{4} = 'issue: add fields genes and grRules to the model';
    checkList{5} = 'not tested: presence of field rules';
    tagReady = 0;
end

% this trick makes the milp problems converge faster
fprintf('4: reducing bigM to improve performance of MILP\n');
indfu = getAllVar(model,{'FU'});
indbu = getAllVar(model,{'BU'});
induf = getAllCons(model,{'UF'});
indur = getAllCons(model,{'UR'});
if full(model.A(induf(1),indfu(1))) < -51
    if length(induf)==length(indur)
        for i = 1:length(induf)
            model.A(induf(i),indfu(i)) = -50;
            model.A(indur(i),indbu(i)) = -50;
        end
    end
    checkList{6} = 'corrected: reduced capacity constrain of fluxes in TFA problem to accelerate MILP convergence';
else
    checkList{6} = 'ok: capacity constrain of fluxes in TFA problem is low and might have been predefined like this by the user';
end

% summarize final status of the model
if isequal(model,modeli) && tagReady
    fprintf('all checks passed: the output model is the same as the input model. The model is ready for phenomapping\n');
elseif ~isequal(model,modeli) && tagReady
    fprintf('checks were applied: the output model is now ready for phenomapping\n');
else
    fprintf('some checks failed: the output model needs manual curation for phenomapping - details in checkList\n');
end

% create temporary results folder to save intermediate and final results: 
% note that PhenoMapping will not upload the intermediate results 
% - but you might need to use them if matlab crashes
if ~isdir(pathToSave(1:end-1))
    mkdir(pathToSave(1:end-1))
end