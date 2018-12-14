function model = addRxnCompSymbol(model)
% Adds the field rxnCompSymbol to a model
%
% USAGE:
%
%    model = addRxnCompSymbol(model)
%
% INPUT:
%    model:           FBA/TFA model structure
%
% OUTPUTS:
%    model:           Model with field rxnCompSymbol
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% check if the compartments for the rxns are defined as in the RAVEN
% toolbox: the conversion is then easy
if isfield(model,'comps') && isfield(model,'rxnComps') && ...
        iscell(model.comps) && isfloat(model.rxnComps)
    model.rxnCompSymbol = model.comps(model.rxnComps);
else
    % we consider the following descriptions to assign compartments to rxns:
    % _compartmentID or [compartmentID] or (compartmentID) at the end of the
    % rxns field in the model
    symbols = {'_','[','('};
    
    % get three last characters of model.rxns field
    rxnCompTag = cell(length(model.rxns),1);
    for i = 1:length(model.rxns)
        rxnCompTag{i} = model.rxns{i}(end-2:end);
    end
    
    % check what tags are part of the three last characters
    numTag = zeros(length(model.rxns),length(symbols));
    for i = 1:length(symbols)
        numTag(:,i) = ~cellfun(@isempty,regexpi(rxnCompTag,symbols{i},'match'));
    end
    
    % assign compartment to rxns based on tag
    rxnCompSymbol = cell(length(model.rxns),1);
    for i = 1:length(model.rxns)
        if sum(numTag(i,:)) == 0
            warning('one rxn does not have a compartment tag. It was assigned to the cytosol. Please check the rxn, the compartment tags of this function and verify this assumption is correct')
            fprintf(strcat('Problematic rxn ID:',model.rxns{i},'\n'));
            rxnCompSymbol{i} = 'c';
        elseif sum(numTag(i,:)) > 1.5
            warning('one rxn has more than one compartment tag at the end. WTF? We assinged the last letter of the string as the compartment. Please check the rxn and the compartment tags of this function and verify this assumption is correct')
            fprintf(strcat('Problematic rxn ID:',model.rxns{i},'\n'));
            rxnCompSymbol{i} = rxnCompTag{i}(end);
        elseif sum(numTag(i,:)) < 1.5
            if strcmp(symbols(numTag(i,:)>0.5),'_')
                rxnCompSymbol{i} = rxnCompTag{i}(end);
            else
                rxnCompSymbol{i} = rxnCompTag{i}(end-1);
            end
        end
    end
    model.rxnCompSymbol = rxnCompSymbol;
end
