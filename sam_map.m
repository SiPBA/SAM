%
% Generate SAM map
%
function [map, p, sigReg] = sam_map(acc, bound, alpha, atlas)

    % Compute z statistic
    [~,p] = sam_zStatistic(acc, bound, alpha);

    % Compute indices of significative regions
    sigReg = p < alpha;
    
    % Initialize SAM map y p-Value map
    map = zeros(size(atlas.nii.img));

    % Estimate all ROIs. Non-releavant ROIs are set to NaN
    for reg = 1:atlas.numReg
        if sigReg(reg), activ = acc(reg); else, activ = nan; end
        map(atlas.nii.img == reg) = activ;
    end
end
