function [union, common, notCommon, num] = getOverlapSets(data1, data2)
% Extract common, not common, union sets for two data sets of same length
% (analysis done at each row)
%
% USAGE:
%
%       [common, notCommon, union, num] = getOverlapSets(data1, data2)
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
    data2 = {};
end

union = data1;
num = zeros(size(data1,1),1);
if ~isempty(data2)
    for i = 1:size(union,1)
        union{i} = unique([data1{i}; data2{i}]);
        num(i) = size(union{i},1);
    end
end

common = union{1};
for i = 2:size(union,1)
    common = common(ismember(common,union{i}));
end

nc = cell(size(union,1),1);
notCommon = [];
for i = 1:size(union,1)
    tmp = union{i,1};
    nc{i} = tmp(~ismember(tmp,common));
    notCommon = unique([notCommon; nc{i}]);
end