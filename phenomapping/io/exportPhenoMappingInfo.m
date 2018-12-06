function exportPhenoMappingInfo(essGenes,conditionToGene,conditionDesc,pathToPhenoData,phenDesc,filename)
% Exports output of PhenoMapping (linkEssGeneIMM2Subs.m, getBotNeckMets)
% to cvs and links genes essential at the conditions with available 
% phenotypes
% 
%
% USAGE:
%
%    exportPhenoMappingInfo(essGenes,conditionToGene,conditionDesc,pathToPhenoData,phenDesc,filename)
%
% INPUT:
%    essGenes:        List of genes that are essential in a condition and  
%                     not in the reference model
%    conditionToGene: Conditions underlying essential genes based on 
%                     phenomapping
%
% OPTIONAL INPUTS:
%    conditionDesc:   Cell array with description of condition (default
%                      = 'condition in PhenoMapping')
%    pathToPhenoData: Cell array with path(s) to phenotype(s) data for the 
%                     organism of study (default = empty / provide a path). 
%                     NOTE: the phenotypes should be saved in a cell called
%                     "phenotypes": 1st and 2nd column geneIDs and 
%                     phenotypes, respectively.
%    phenDesc:        Cell array with description of phenotype(s) (default
%                      = 'phenotype')
%    filename:        Name used to save file (default = 
%                     'PhenoMappingCondition2Genes'). NOTE: no need to
%                     provide the extension of the file here. It will be
%                     csv by default. But the path where the file will be
%                     saved should be provided else it will be saved in the
%                     folder where you currently are.
%
% OUTPUTS:
%    cvs:             Saved in phenomapping/tmpresults. You can open the  
%                     cvs in excel and separate Text to columns with Tab
%                     as a delimiter.
%
% .. Author:
%       - Anush Chiappino-Pepe 2018
%

if (nargin < 3) || isempty(conditionDesc)
    conditionDesc = {'condition in PhenoMapping'};
end
if (nargin < 4)
    pathToPhenoData = [];
end
if (nargin < 5) || isempty(phenDesc)
    phenDesc = pathToPhenoData;
end
if (nargin < 6) || isempty(filename)
    filename = 'PhenoMappingCondition2Genes';
end

if (size(phenDesc,1)>size(phenDesc,2))
    phenDesc = phenDesc';
end

[conditionToGene] = collectInfoCell(conditionToGene);

NumAlt = size(conditionToGene{1,1},2);
heading = cell(1,NumAlt);
for i = 1:NumAlt
    heading(i) = strcat('Alt',{' '},num2str(i));
end

data = [];
for i = 1:size(conditionToGene,1)
    if ~isempty(conditionToGene{i})
        conditionToGene{i} = strrep(conditionToGene{i},strcat({' '},'<=>'),'');
        conditionToGene{i} = strrep(conditionToGene{i},',','');
        if ~isempty(pathToPhenoData)
            phen = [];
            for j = 1:length(pathToPhenoData)
                tmp = getPhenotype(essGenes{i},'','',pathToPhenoData{j});
                phen = [phen, tmp(:,2)];
            end
        else
            phenDesc = {'phenotype'};
            tmp = getPhenotype(essGenes{i},'','',[]);
            phen = tmp(:,2);
        end
        data = [data; ...
            [strcat(conditionDesc,{' '},num2str(i)),phenDesc,heading]; ...
            [essGenes{i},phen,conditionToGene{i}]];
    end
end

desc = '%s\t';
for i = 2:size(data,2)
    desc = strcat(desc,'%s\t');
end

writeData(strcat(filename,'.csv'),data,desc);