function [model, checkList, tagReady] = initTestPhenoMappingModel(modeli)
% Initial test of the model to check if it is ready for a PhenoMapping
% analysis
%
% USAGE:
%
%       [model, checkList, tagReady] = initTestPhenoMappingModel(modeli)
%
% INPUTS:
%    model:           model with FBA/TFA structure
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

tagReady = 1;
model = modeli;

fprintf('1: checking model structure\n');
if isfield(model,'A') && isfield(model,'f') && isfield(model,'var_lb') ...
        && isfield(model,'var_ub') && isfield(model,'rhs') ...
        && isfield(model,'constraintNames') && isfield(model,'varNames') ...
        && isfield(model,'constraintType') && isfield(model,'vartypes') ...
        && isfield(model,'objtype')
    checkList{1} = 'ok: the model already has a TFA structure';
else
    model = convGEMToTFAStruc(model);
    checkList{1} = 'corrected: the model was converted to a TFA structure';
end

fprintf('2: checking fields for gene essentiality\n');
if isfield(model,'genes') && isfield(model,'grRules')
    checkList{2} = 'the model has the fields for gene essentiality';
    if isfield(model,'rules')
        checkList{3} = 'ok: the field rules is present';
    else
        [model] = generateRules(model);
        checkList{3} = 'corrected: the field rules was added';
    end
else
    checkList{2} = 'issue: add fields genes and grRules to the model';
    checkList{3} = 'not tested: presence of field rules';
    tagReady = 0;
end

indfu = getAllVar(model,{'FU'});
indbu = getAllVar(model,{'BU'});
induf = getAllCons(model,{'UF'});
indur = getAllCons(model,{'UR'});
if model.A(induf(1),indfu(1)) > -51
    if length(induf)==length(indur)
        for i = 1:length(induf)
            model.A(induf(i),indfu(i)) = -50;
            model.A(indur(i),indbu(i)) = -50;
        end
    end
    checkList{4} = 'corrected: reduced capacity constrain of fluxes in TFA problem to accelerate MILP convergence';
else
    checkList{4} = 'ok: capacity constrain of fluxes in TFA problem is low and might have been predefined like this by the user';
end

if isequal(model,modeli) && tagReady
    fprintf('all checks passed: the output model is the same as the input model. The model is ready for phenomapping\n');
elseif ~isequal(model,modeli) && tagReady
    fprintf('checks were applied: the output model is now ready for phenomapping\n');
else
    fprintf('some checks failed: the output model needs manual curation for phenomapping - details in checkList\n');
end

end