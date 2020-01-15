% Extract a summaryzed set of feature from a dataset
%
% Usage:
%   feats = sam_featureExtraction(data, group, method, nFeats)
% Inputs:
%   data:      Data matrix (Rows: observations, Columns: variables).
%   group:     Labels for observations in 'data'.
%   method:    Posible values: 'pls' for Partial Least Squares.
%   nFeat:     Number of features extracted.
% Outputs:
%   feats:     Feature matrix (Rows: observations, Columns: features).
%
% Please cite:
%   Juan M. Gorriz et al., A Machine Learning Approach to Reveal the 
%   Neuro-Phenotypes of Autisms, International Journal of Neural Systems, 
%   doi: 10.1142/S0129065718500582
%
function feats = sam_featureExtraction(data, group, method, nFeats)
    switch (method)
        case 'pls'
            [~,~,feats] = plsregress(zscore(data), group, nFeats);
        otherwise
            error([method ': Method not implemented yet']);
    end
end
