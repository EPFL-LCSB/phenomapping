function [levelGenes, genesNA] = generateGeneLevels(model, genesrnainfo,...
    genesrnainfovalue, tagEstNA)
% Generate the gene levels field input to texfba for analysis of
% genome-scale models with transcriptomics data
%
% [levelGenes, genesNA] = generateGeneLevels(model, genesrnainfo, genesrnainfovalue, tagEstNA)
%
% INPUT
% model             FBA model structure
% genesrnainfo      List of genes with associated mRNA values
% genesrnainfovalue List of mRNA cound values for the genes in genesrnainfo
%                   - if a matrix is provided the mean of all values for
%                   the same gene will be calculated
%
% OPTIONAL INPUTS:
% tagEstNA:         True to estimate for genes without data the average
%                   value of mRNA level - if the hl and hp 
%                   parameters of the texfba methodology are far from the 
%                   mean this means these genes will be neglected anyways. 
%                   False to define them and nan (default = false) 
%
% OUTPUTS
% levelGenes        Levels of mRNA associated to the genes in the model (in
%                   the order as they appear in model.genes)
% genesNA           Genes in the model for which no data was found in the
%                   dataset and values assigned based on tagEstNA
%
%
% .. Author:
% Anush Chiappino-Pepe 2016
% 

if (nargin < 4)
    tagEstNA = 0;
end

levelGenes = zeros(length(model.genes),1);

% preprocessing of genes in model and genesrnainfo
genesrnainfo = strtok(genesrnainfo,'.');
genesrnainfo = strrep(genesrnainfo,{' '},'');
genes = strrep(model.genes,'putative_','');
genes = strtok(genes,'.');
genes = strrep(genes,{' '},'');
    
% find genes in the model in the input genelist with transcriptomics data 
[y,row] = ismember(genes, genesrnainfo);

levelGenes(y) = mean(genesrnainfovalue(row(y),:),2);
if tagEstNA
    levelGenes(~y) = mean(mean(levelGenes(y),1),2); % for no data assume the mean of all - this will eliminate the gene from the analysis
else
    levelGenes(~y) = nan;
end

% extract info of genes without data
if any(~y)
    genesNA = [model.genes(~y), genes(~y), levelGenes(~y)];
else
    genesNA = {'data for all genes was found'};
end