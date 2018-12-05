function [model,deletedRxns] = thermoDeleteModelRxns(model,rxnList,downRegFraction)
%Delete one or more rxns and appends '_deleted' to the rxn(s)
%
% [model,deletedRxns] = thermoDeleteModelGenes(model,rxnList,downRegFraction)
%
%INPUT
% model             TFA model structure
%
%OPTIONA INPUTS
% rxnList           List of rxns to be deleted (Default = empty)
% downRegFraction   Fraction of the original bounds that the reactions
%                   will be assigned
%                   (Default = 0 corresponding to a full deletion)
%
%OUTPUTS
% model             TFA model with the selected genes deleted
% deletedRxns       The list of rxns removed from the model.  
%

if (nargin < 2)
    rxnList = [];
end

if (nargin < 3)
    downRegFraction = 0;
end

deletedRxns = rxnList(ismember(rxnList,model.rxns));

if (nargin > 2)
    model = changeRxnBounds(model,rxnList,downRegFraction*model.lb(findRxnIDs(model,rxnList)),'l');
    model = changeRxnBounds(model,rxnList,downRegFraction*model.ub(findRxnIDs(model,rxnList)),'u');
    model.var_ub(ismember(model.varNames,strcat('R_', rxnList)))=downRegFraction*model.var_ub(ismember(model.varNames,strcat('R_', rxnList)));
    model.var_ub(ismember(model.varNames,strcat('F_', rxnList)))=downRegFraction*model.var_ub(ismember(model.varNames,strcat('F_', rxnList)));
    model.var_lb(ismember(model.varNames,strcat('R_', rxnList)))=downRegFraction*model.var_lb(ismember(model.varNames,strcat('R_', rxnList)));
    model.var_lb(ismember(model.varNames,strcat('F_', rxnList)))=downRegFraction*model.var_lb(ismember(model.varNames,strcat('F_', rxnList)));
    model.var_lb(ismember(model.varNames,strcat('NF_', rxnList)))=downRegFraction*model.var_lb(ismember(model.varNames,strcat('NF_', rxnList)));
    model.var_ub(ismember(model.varNames,strcat('NF_', rxnList)))=downRegFraction*model.var_ub(ismember(model.varNames,strcat('NF_', rxnList)));
else
    %                 Full deletion
    model = changeRxnBounds(model,rxnList,0,'b');
    model.var_ub(ismember(model.varNames,strcat('R_', rxnList)))=0;
    model.var_ub(ismember(model.varNames,strcat('F_', rxnList)))=0;
    model.var_lb(ismember(model.varNames,strcat('R_', rxnList)))=0;
    model.var_lb(ismember(model.varNames,strcat('F_', rxnList)))=0;
    model.var_lb(ismember(model.varNames,strcat('NF_', rxnList)))=0;
    model.var_ub(ismember(model.varNames,strcat('NF_', rxnList)))=0;
end

