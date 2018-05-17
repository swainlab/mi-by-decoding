function I = Info(M)
%INFO Information estimation from a confusion matrix
%   Info(M) returns the estimated mutual information (in bits) for NxN 
%   confusion matrix M
%
%   If you publish results that make use this software or the mutual 
%   information by decoding algorithm, please cite:
%   Granados, A.A., Pietsch, J.M.J., Cepeda-Humerez, S.A., Farquhar, I.L.,
%   Tkacik, G., and Swain, P.S. (2018) Distributed and dynamic
%   intracellular organization of extracellular information.
%
%   Authors: Sarah Cepeda, Julian Pietsch
%   Copyright Gasper Tkacik and Peter Swain 2018

if any(M(:)<0)
    error('MutualInformationCode:Info:NegativeConfMatrix',...
        'The confusion matrix cannot contain negative values');
end
if any(M(:)==0)
    % Add a small eps for numerical stability when probabilities approach 
    % zero (valid approximation by L'Hopital's rule for the calculation of 
    % mutual information below):
    M = M + min(M(M(:)~=0))*10^-8;
end

[len_M,~] = size(M);
L  = log2(M);
P1 = sum(M,1);
P2 = sum(M,2);

H = M.*L; H = H(~isnan(H));
% I = sum_x sum_y P(x,y) * log2(P(x,y)/(P(x)*P(y)))
I = sum(H(:)) - sum(sum(M.*log2(repmat(P1,len_M,1).*repmat(P2,1,len_M))));

if sum(M(:)) < 0.00009
    error('MutualInformationCode:Info:NegativeConfMatrix', ...
        'The confusion matrix contains negative values');
end

end
