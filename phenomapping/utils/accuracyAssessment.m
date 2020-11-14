function [MCC, accuracy, NPR, PPR, sensitivity, specificity] = ...
    accuracyAssessment(geneList, inSilicoPheno, inVivoPheno)
% Computes the accuracy metric for a list of genes with mapped in silico
% and in vivo phenotypes
%
% USAGE:
%
%    [MCC, accuracy, NPR, PPR, sensitivity, specificity] = accuracyAssessment(geneList, inSilicoPheno, inVivoPheno)
%
% INPUT:
%    geneList:        Gene list
%    inSilicoPheno:   In silico phenotypes for the gene list as obtained
%                     with an essentiality prediction. Essential should be
%                     defined as 'essential' and non-essential as
%                     'dispensable'.
%    inVivoPheno:     In vivo phenotypes for the gene list as obtained
%                     experimentally. Essential should be defined as 
%                     'essential' and non-essential as 'dispensable'.
%                     
%
% OUTPUTS:
%    MCC:             Matthew Correlation Coefficient
%    accuracy:        Overall Accuracy
%    NPR:             Negative Prediction Rate
%    PPR:             Positive Prediction Rate
%    sensitivity:     Sensitivity
%    specificity:     Specificity
%
% .. Author:
% Anush Chiappino-Pepe 2020
%

% lower case
inSilicoPheno = lower(inSilicoPheno);
inVivoPheno = lower(inVivoPheno);

positiveTag = 'dispensable';
negativeTag = 'essential';

if ~any(ismember(positiveTag,inSilicoPheno))
    warning('there is no dispensable phenotype in silico')
end
if ~any(ismember(positiveTag,inVivoPheno))
    warning('there is no dispensable phenotype in vivo')
end
if ~any(ismember(negativeTag,inSilicoPheno))
    warning('there is no essential phenotype in silico')
end
if ~any(ismember(negativeTag,inVivoPheno))
    warning('there is no essential phenotype in vivo')
end

% extract negatives and positives
inSilicoNeg = geneList(ismember(inSilicoPheno, negativeTag));
inVivoNeg = geneList(ismember(inVivoPheno, negativeTag));

inSilicoPos = geneList(ismember(inSilicoPheno, positiveTag));
inVivoPos = geneList(ismember(inVivoPheno, positiveTag));

% identify TPs, TNs, FPs, and FNs
tn = sum(ismember(inSilicoNeg, inVivoNeg));
tp = sum(ismember(inSilicoPos, inVivoPos));
fp = sum(ismember(inSilicoPos, inVivoNeg));
fn = sum(ismember(inSilicoNeg, inVivoPos));

% Metrics for accuracy assessment
% Matthew correlation coefficient
MCC = (tp*tn - fp*fn) / (sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)));

% overal accuracy
accuracy = (tp + tn) / (tp + fp + tn +fn);

% negative prediction rate
NPR = tn / (tn + fn);

% positive prediction rate
PPR = tp / (tp + fp);

% sensitivity
sensitivity = tp / (tp + fn);

% specificity
specificity = tn / (tn + fp);

end

