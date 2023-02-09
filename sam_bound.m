%
% Calculate upper bound for the classification error.
%
function res = sam_bound(varargin)
    f = eval(['@bound_' varargin{1}]);
    res = f(varargin{2:end});
end

function l = bound_getList(varargin)
    l = {'G','GZ','Vapnik','Bayes'};
end


% Method 1: Vapnik´s upper bound
function bound = bound_Vapnik(n, d, alpha, varargin)
    bound = sqrt(abs(((d+1)*(log(2*n/(d+1))+1)-log(alpha/4))./n));
end



% Method 2: igp upper bound
function bound = bound_G(n, d, alpha, varargin)
    cld=0;
    for k=1:d, cld = cld + 2 * nchoosek(n-1,k-1); end
    bound = sqrt(log(cld/alpha)/(2*n));

end



% Method 3: igp upper bound
function bound = bound_GZ(n, d, alpha, varargin)
    cld=0;
    for z=0:(d-1)
        for k=1:d-z
            cld = cld + 2*nchoosek(n,z)*nchoosek(n-1-z,k-1);
        end
    end
    bound = sqrt(log(cld/alpha)/(2*n));
end


function bound = bound_Bayes(n, d, alpha, Mdl, Acc_emp, dropout, varargin)
    lambda=1/2+0.1:0.01:10; % parameter upper bounding the root
    a=1./(1-1./(2*lambda));
    
    % For each standardized region one bound...
    Lmax=1; % maximum error
    %Bound: L=Lemp+Bound is the worst case, Bound=(a-1)*Lemp+a*B
    %where B= (Lambda*Lmax)/n * ( (1-dropout)/2 ||theta||^2 + ln (1/alpha) )
    for i=1:numel(Acc_emp)
        
        Theta=[Mdl{i}.BinaryLearners{1}.Beta Mdl{i}.BinaryLearners{1}.Bias]; % vector of parameters
        
        bound(i)=min( (a-1)*(1-Acc_emp(i))+...
            a.*(lambda*Lmax/n)*( (1-dropout)/2 * ...
            norm(Theta)^2 + log(numel(lambda)/alpha) ) ); % minimum in lambda
    end
end

