% Tests for the phenomapping core transcriptomics module
%
% USAGE:
%
%    run('test_core_transcriptomics.m')
%
%
% .. Author:
% Anush Chiappino-Pepe 2018
%

% run test_core.m before or uncomment this part
% clear
% clc
% 
% this_file_directory = mfilename('fullpath');
% run(strrep(this_file_directory,'test_core_transcriptomics','test_core.m'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENE ESSENTIALITY WITH TFA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description: Get essentiality of the model in the reference conditions  
% (normally in rich medium and without any data integrated)

% inputs
grRate = 0.05;          % optimal growth rate
essThr = 0.1;           % essentiality threshold

% analysis
if ~exist('essTFAref','var') || isempty(essTFAref)
    [~, grRateGeneTFA] = thermoSingleGeneDeletion(model, 'TFA', ...
        model.genes, 0, 0, 0, essThr, model.indNF);
    grRateGeneTFA(isnan(grRateGeneTFA)) = 0; %by default NaN is considered an essential KO
    essTFAref = model.genes(grRateGeneTFA < essThr*grRate);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRANSCRIPTOMICS ANALYSIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs




save(strcat(saving_directory,filename,'.mat'));

