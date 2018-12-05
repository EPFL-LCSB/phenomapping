function exportIMM2ess2subInfo(essIMMaddToIRM,subsToGenes,filename,pathToData,phenotypeDesc)
% Exports output of linkEssGeneIMM2Subs.m to cvs and links genes essential 
% at the IMM with available phenotypes
%
% USAGE:
%
%    exportIMM2ess2subInfo(essIMMaddToIRM,subsToGenes,filename,pathToData,phenotypeDesc)
%
% INPUT:
%    essIMMaddToIRM:  Essential genes in IMM on the top of essential genes
%                     in rich medium
%    subsToGenes:     Substrates whose abscence in the medium determine
%                     the essential function of the genes in essIMMaddToIRM
%
% OPTIONAL INPUTS:
%    filename:        Name used to save file (default = 'imm_ess2sub')
%    pathToData:      Cell array with path(s) to phenotype(s) for a life 
%                     stage of the organism of study (default = empty / get
%                     data for P. berghei blood and liver stages saved in 
%                     phenomapping/data/pbe). NOTE: new phenotypes should
%                     be saved in a cell called "phenotypes".
%    phenotypeDesc:   Cell array with description of phenotype(s) (default
%                      = 'phenotype')
%
% OUTPUTS:
%    cvs:             Saved in phenomapping/tmpresults
%
% .. Author:
%       - Anush Chiappino-Pepe 2018

% Ways to call this function
% provide only the outputs of linkEssGeneIMM2Subs.m
% exportIMM2ess2subInfo(essIMMaddToIRM,subsToGenes);
% provide a file name
% exportIMM2ess2subInfo(essIMMaddToIRM,subsToGenes,strcat(filename,'_ess2subjoint'));
% provide path to one or more sets of phenotypes to map to the genes
% essential at the imm
% exportIMM2ess2subInfo(essIMMaddToIRM,subsToGenes,[],{'/phenomapping/data/pbe/pbe_phenotypes_blood_Nov16.mat';'/phenomapping/data/pbe/pbe_phenotypes_liver_Jul18.mat'},{'blood pbe';'liver pbe'});
% 

if (nargin < 3) || isempty(filename)
    filename = 'imm_ess2sub';
end
if (nargin < 4)
    pathToData = [];
end
if (nargin < 5)
    phenotypeDesc = cell(1,length(pathToData));
    phenotypeDesc(:) = {'phenotype'};
end

if (size(phenotypeDesc,1)>size(phenotypeDesc,2))
    phenotypeDesc = phenotypeDesc';
end

[subsToGenes] = collectInfoCell(subsToGenes);

NumAlt = size(subsToGenes{1,1},2);
heading = cell(1,NumAlt);
for i = 1:NumAlt
    heading(i) = strcat('Alt',{' '},num2str(i));
end

data = [];
for i = 1:size(subsToGenes,1)
    if ~isempty(subsToGenes{i})
        subsToGenes{i} = strrep(subsToGenes{i},strcat({' '},'<=>'),'');
        subsToGenes{i} = strrep(subsToGenes{i},',','');
        if isempty(pathToData)
            phenotypeDesc = {'blood phenotype','liver phenotype'};
            b = getPhenotype(essIMMaddToIRM{i},'blood','pbe');
            l = getPhenotype(essIMMaddToIRM{i},'liver','pbe');
            phen = [b(:,2), l(:,2)];
        else
            phen = [];
            for j = 1:length(pathToData)
                tmp = getPhenotype(essIMMaddToIRM{i},'','',pathToData{j});
                phen = [phen, tmp(:,2)];
            end
        end
        data = [data; ...
            [strcat('IMM',{' '},num2str(i)),phenotypeDesc,heading]; ...
            [essIMMaddToIRM{i},phen,subsToGenes{i}]];
    end
end

desc = '%s\t';
for i = 2:size(data,2)
    desc = strcat(desc,'%s\t');
end

writeData(strcat('tmpresults/',filename,'.csv'),data,desc);