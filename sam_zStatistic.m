% Computing Z-statistic, p-values and threshold
%
% Usage:
%   [thres, p, zAcc, interval] = z_statistic(Acc, Pi0, alpha)
% Input:
%   acc:      Accuracy values (proportions)
%   bound:    Analytical upper bound 
%   alpha:    Significance level
% Outputs:
%   thres:    Threshold for significance
%   p:        p value
%   zAcc:     Z-Statistic for a proportion test 
%   interval: Confidence interval 
%
% Please cite:
%   J.M. Gorriz et al. for the Alzheimer's Disease Neuroimaging Initiative 
%   (ADNI) and for the Parkinson's Progression Markers Initiative (PPMI). 
%   STATISTICAL AGNOSTIC MAPPING: A FRAMEWORK IN NEUROIMAGING BASED ON 
%   CONCENTRATION INEQUALITIES. 
%   doi: https://doi.org/10.1101/2019.12.27.889436.
%
function [thres, p, zAcc, interval] = sam_zStatistic(acc, bound, alpha)

    accCor = (acc - bound)';       % Accuracy values (corrected)
    pi0 = mean([accCor; acc']);  % Mean value for the known proportions 
                                 % (worst case corrected and uncorrected) 

    % Accuracies in columns
    ind = accCor > 0.5;  % Accuracies better than random, i.e. pi>0.5
    %ind = accCor > mean(acc);  % Accuracies better than the mean.
    l = sum(ind);     % It considers only accuracies in confidence interval given by concentration inequalities
    %l = numel(Accn); % It considers all the accuracies, less conservative

    se0 = sqrt(pi0*(1-pi0)/l);   % standard error for a proportion
    zAcc = (accCor-pi0)/se0;     % z-statistic for a proportion

    % Normality assumed on H0
    p = 1 - normcdf(zAcc,mean(zAcc(ind)),std(zAcc(ind)),'alpha',alpha); 

    % Invert PDF assuming Z_Acc gaussian (risky)
    thres = norminv(p,mean(zAcc),std(zAcc));                           

    % Computing interval of this test
    interval=[accCor-1.96*sqrt(accCor.*(1-accCor)/l) 
              accCor+1.96*sqrt(accCor.*(1-accCor)/l)];

end
