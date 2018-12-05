function model = loadConstraints(model,constraints,type)
% loads the constraints into the model
% the constraint structure MUST have the following fields:
% first column specifies the variable names as per the varNames
% second column specifies the lower bound
% third colum specifies the upper bound

if nargin < 3
    type = 'TFA';
end

if strcmp(type,'TFA')
    [num_vars tmp] = size(constraints);

    for i=1:num_vars
        varName = constraints{i,1};
        var_index = find(ismember(model.varNames,varName));
        if isempty(var_index)
            fprintf('%s not found\n',constraints{i,1});
        else
            fprintf('%s lb: %d => %d\tub: %d => %d\n',constraints{i,1},model.var_lb(var_index),constraints{i,2},model.var_ub(var_index),constraints{i,3}); 

            model.var_lb(var_index) = constraints{i,2};
            model.var_ub(var_index) = constraints{i,3};
        end
    end
elseif strcmp(type,'FBA')
    [num_vars tmp] = size(constraints);
    
    for i=1:num_vars
        varName = constraints{i,1};
        var_index = find(ismember(model.rxns,varName));
        
        if isempty(var_index)
            fprintf('%s not found\n',constraints{i,1});
        else
            fprintf('%s lb: %d => %d\tub: %d => %d\n',constraints{i,1},model.lb(var_index),constraints{i,2},model.ub(var_index),constraints{i,3}); 

            model.lb(var_index) = constraints{i,2};
            model.ub(var_index) = constraints{i,3};
        end
    end
end