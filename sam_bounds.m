% Computing analytical upper bounds
%
% Usage:
%   [boundVC, boundGZ, boundG] = sam_bounds(l, d, alpha)
% Inputs:
%   n:       Sample size
%   d:       Dimensions
%   alpha:   Significance level
% Outputs:
%   boundVC: Vapnik bound
%   boundG:  Method 1 upper bound
%   boundGZ: Method 2 upper bound
%
% Please cite:
%   Gorriz, J.M. et al. On the computation of distribution-free performance 
%   bounds: Application to small sample sizes in neuroimaging. 
%   Pattern Recognition 93, 1-13
% Method 1:
%   V. Vapnik Estimation dependencies based on Empirical Data.
%   Springer-Verlach. 1982 ISBN 0-387-90733-5
% Method 2&3:
%   Górriz, J.M. et al. On the computation of distribution-free performance bounds: 
%   Application to small sample sizes in neuroimaging. 
%   Pattern Recognition 93, 1-13
%
function [boundVC, boundG, boundGZ] = sam_bounds(n, d, alpha)

    % Method 1: Vapnik´s upper bound
    boundVC = sqrt(abs(((d+1)*(log(2*n/(d+1))+1)-log(alpha/4))./n));


    % Method 2: igp upper bound
    cld=0;
    for k=1:d, cld = cld + 2 * nchoosek(n-1,k-1); end
    boundG = sqrt(log(cld/alpha)/(2*n));

    % Method 3: igp upper bound
    cld=0;
    for z=0:(d-1)
        for k=1:d-z
            cld = cld + 2*nchoosek(n,z)*nchoosek(n-1-z,k-1);
        end
    end
    boundGZ = sqrt(log(cld/alpha)/(2*n));

end


