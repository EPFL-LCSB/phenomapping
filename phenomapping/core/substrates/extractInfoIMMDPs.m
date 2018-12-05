function [StatsMets, Mets, drainClass, statsClass] = ...
    extractInfoIMMDPs(model, DPs, indUSE)
% Extract statistics of metabolites in DPs of IMM/IMS
%
% USAGE:
%
%    [StatsMets, Mets, drainClass, statsClass] = extractInfoIMMDPs(model, DPs, indUSE)
%
% INPUT:
%    model:           TFA model with MILP structure for IMM analysis
%    DPs:             Directionality profile matrix with alternatives in
%                     each column from IMM analysis
%
% OPTIONAL INPUTS:
%    indUSE:          indexes of integers used in the MILP (default = model.indUSE)
%
%
% OUTPUTS:
%    StatsMets:       Statistics on appearence of mets among all
%                     alternatives
%    Mets:            Summary of drains and metabolites
%    drainClass:      Elements studied and their classification
%    statsClass:      Statistics on classification
%
% .. Author:
% Anush Chiappino-Pepe 2014
% 

% Extracting info only for integer variables in DPs
if nargin < 3 || isempty(indUSE)
    indUSE = model.indUSE;
end
InfoDPs = DPs(indUSE,:);
drains = model.varNames(indUSE);

% change DP matrix to be 1 if element is active
InfoDPs = ones(size(InfoDPs,1),size(InfoDPs,2))-InfoDPs;

NumActive = zeros(1,size(DPs,2));
% identify number of mets in the IMM of each alternative
for i = 1:size(InfoDPs,2)
    NumActive(i) = round(sum(InfoDPs(:,i)),2);
end
if ~isequal(min(NumActive),max(NumActive))
    InfoDPs = InfoDPs(:,NumActive<min(NumActive)+0.5);
    fprintf('check: not all DPs have same IMM size\n');
end

% classify the integer variables
StatsMets = zeros(length(drains),1);
drainClass = cell(size(InfoDPs,1),1);
statsClass = {'constitutive'; 'not part of alternatives'; 'non-constitutive'; 'size of alternatives'; 'min number of groups'; 'number of alternatives'};
numClass = zeros(6,1);

% get statistics on appearence of mets among all alternatives
for i = 1:size(InfoDPs,1)
    StatsMets(i) = round(sum(InfoDPs(i,:)),2);
    if (sum(InfoDPs(i,:)) > (size(InfoDPs,2)-0.5))
        drainClass(i,1) = {'constitutive'};
        numClass(1) = numClass(1) + 1;
    elseif (sum(InfoDPs(i,:)) < 1)
        drainClass(i,1) = {'not part of alternatives'};
        numClass(2) = numClass(2) + 1;
    else
        drainClass(i,1) = {'non-constitutive'};
        numClass(3) = numClass(3) + 1;
    end
end
numClass(4) = sum(InfoDPs(:,1));
numClass(5) = numClass(4) - numClass(1);
numClass(6) = size(InfoDPs,2);
statsClass(:,2) = num2cell(numClass);

% get met names associated to statistics (should be also equivalent to the
% output drains of analysisIMM.m)
drains = strrep(drains,'BFUSE_','');
drains = strrep(drains,'R_','');
drains = strrep(drains,'F_','');
Mets = drains;
Mets(:,2) = printRxnFormula(model,drains,0);
Mets(:,3) = printRxnFormula(model,drains,0,0,1);

end