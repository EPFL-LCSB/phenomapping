function [DPs, model, objectives] = findAltConc(model, NumAlt, indUSE, time)
% Get alternative solutions for MILP (minimization) concentrations
%
% USAGE:
%
%       [DPs, model, objectives] = findAltConc(model, NumAlt, indUSE, time)
%
% INPUTS:
%    model:           Model with TFA structure and MILP formulation
%
% OPTIONAL INPUTS:
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    intTag:          Integet tag used in the MILP (default = 'BFUSE')
%    time:            Time in sec after which simulation will be stopped
%                     (default = empty / no time provided)
%
% OUTPUTS:
%    DPs:             Directionality profile matrix with alternatives in
%                     each column
%    model:           Model with integer cuts integrated to avoid
%                     repetition of same solution or supersolution (this 
%                     is a superset of an obtained solution)
%    objectives:      Optimal value in MILP, here number of LC constraints
%                     that couldnot be integrated into the model
%
% .. Author:
% Meric Ataman & Anush Chiappino-Pepe 2015
% 
if (nargin < 2)
    NumAlt = 1;
end
if (nargin < 3) || isempty(indUSE)
    intTag = {'LCUSE'};
    indUSE = getAllVar(model,intTag);
end
if (nargin < 4)
    time = [];
end

[~,NumVars] = size(model.A);
model.f = zeros(NumVars,1);
model.f(indUSE) = 1;
model.objtype = 1;

NumSols = 0;
sol = optimizeThermoModel(model,0,'cplex',time);
DPs = [];
objectives = [];

while ((NumSols < NumAlt) && ~(isempty(sol.x)))
    
    [NumCons,NumVars] = size(model.A);
    
    if ~(isempty(sol.x))
        NumSols = NumSols + 1;
        objectives(NumSols,1) = sol.val;
        DPs(:,NumSols) = sol.x;
        USEvecSol = zeros(NumVars,1);
        USEvecSol(indUSE) = sol.x(indUSE);
        actUSEvec = find(USEvecSol > 0.98);
        
        cardinality = length(actUSEvec);
        NewCons = zeros(NumVars,1);
        NewCons(actUSEvec) = 1;
        
        model.A(NumCons+1,:) = NewCons;
        model.rhs(NumCons+1) = cardinality - 0.5;
        model.constraintNames{NumCons+1} = ['CUT_' num2str(NumSols)];
        model.constraintType{NumCons+1} = '<';
        
        % solve the problem again
        sol = optimizeThermoModel(model,0,'cplex',time);
        
        fprintf('Number of DPs:\t%d\n',NumSols);
        
        if rem(NumSols,10) == 0
            save DPs DPs;
        end
    end
end
end




