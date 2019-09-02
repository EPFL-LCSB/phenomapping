function allTasksImpact = extractTasksImpact(impactTasks, geneList)
% Analyze output from checkBBBTasks
%
% USAGE:
%
%       [union, common, notCommon, num] = getOverlapSets(data1, data2)
%
% INPUTS:
%    data1:           Data set 1 - row
%
% OPTIONAL INPUTS:
%    data2:           Data set 2 - row (default = empty)
%
% OUTPUTS:
%    common:          data common between data sets 1 and 2
%    notCommon:       data differing between data sets 1 and 2
%    union:           joint data 1 and 2 at each cell (row)
%    num:             length of each cell (row)
%
% .. Author:
% Anush Chiappino 2015
% 

if (nargin < 2)
    geneList = {'infesible model'};
end

allTasksImpact = cell(length(geneList),3);

for i = 1:size(impactTasks,1)
    allTasksImpact(i,1) = geneList(i);
    if ~isempty(impactTasks{i})
        id = impactTasks{i}(1,1);
        name = impactTasks{i}(1,2);
        for j = 2:size(impactTasks{i},1)
            id = strcat(id,'|', impactTasks{i}(j,1));
            name = strcat(name,'|', impactTasks{i}(j,2));
        end
        allTasksImpact(i,2) = id;
        allTasksImpact(i,3) = name;
    else
        allTasksImpact(i,2) = {''};
        allTasksImpact(i,3) = {''};
    end
end