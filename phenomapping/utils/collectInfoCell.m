function [InfoCollected] = collectInfoCell(InfoSubCell)
% Collect info of subcells in one main cell
%
% USAGE:
%
%    InfoCollected = collectsubstrates(InfoSubCell)
%
% INPUT:
%    InfoSubCell:     Cell with data in subcells
%
%
% OUTPUTS:
%    InfoCollected:   Cell with all data collected
%    indCell:         indexes of cells in which 
%
% .. Author:
% Anush Chiappino-Pepe 2017
% 
% NOTE: This might not be the most efficient way to do this, but it
% works well for the output from linkEssGeneIMM2Subs.m
%

% InfoCollected is supposed to be a column matrix with cells
InfoCollected = InfoSubCell;
if (size(InfoCollected,1) < size(InfoCollected,2))
    InfoCollected=InfoCollected';
end

for i = 1:size(InfoCollected,1) % go throw each row of InfoCollected (a cell)
    if (~isempty(InfoCollected{i,1})) % that is not empty
        for j = 1:size(InfoCollected{i,1},1) % go throw each row of this cell
            for z = 1:size(InfoCollected{i,1},2) % go throw each column of this cell
                if (~isempty(InfoCollected{i,1}{j,z})) % that is not empty
                    % concatenate all elements in this cell
                    InfoCollected{i,1}{j,z} = concatenateList(InfoCollected{i,1}{j,z}, '|');
                end
            end
        end
    end
end
end