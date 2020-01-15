% Rank features and select the highest ranked ones
%
% Usage:
%   [featIdx, featRank] = featureSelection(data, group, method, maxFeat)
% Inputs:
%   data:      Data matrix (Rows: observations, Columns: variables).
%   group:     Labels for observations in 'training'.
%   method:    Posible values: 'anova', 'ttest','entropy','wilcoxon','roc'
%              and 'bhattacharyya'.
%   nFeat:     Number of features seleted.
% Outputs:
%   featIdx:   Indices of selected features.
%   featRank:  Rank of selected features.
%
% Please cite:
%   Juan M. Gorriz et al., A Machine Learning Approach to Reveal the 
%   Neuro-Phenotypes of Autisms, International Journal of Neural Systems, 
%   doi: 10.1142/S0129065718500582
%
function [featIdx, featRank] = sam_featureSelection(data, group, method, nFeat)
    if strcmp(method, 'anova')
        for f=1:size(data,2)
            [p, tbl, stats] = anova1(data(:,f)',group','off');
            matF(f)=tbl{2,5};
        end
        [featRank, featIdx] = sort(matF,'descend');
        featIdx = featIdx(~isnan(featRank));
        featRank = featRank(~isnan(featRank));
    else
        [IDX, Z] = rankfeatures(data', group, 'CRITERION', method);
        featIdx = IDX;
        featRank = Z(IDX);
    end
        
    nFeat = min(nFeat, numel(featIdx));
    featIdx = featIdx(1:nFeat);
    featRank = featRank(1:nFeat);
end
