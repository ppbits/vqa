%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporal Pooling Method
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Nov 17, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dv = TemporalPooling(D, lambda1, beta, lambda2)
if nargin < 4
    lambda2 = 0.5;
end
if nargin < 3
    beta = 1; %mean
end
if nargin < 2
    lambda1 = 1.5;
end

dv = StandardDiviationShiftedMean(D, lambda1, beta, lambda2);
end

function dv = StandardDiviationShiftedMean(D, lambda, beta, lambda2)
if nargin < 4
    thres1 = 1.5;
    thres2 = 0.5;
else
    thres1 = 1+lambda2;
    thres2 = 1-lambda2;
end

D =D(:);
deta = std(D);
u = mean(D);
alpha = 1 + lambda*deta/u;

if lambda > 0
    if alpha > thres1
        alpha = thres1;
    end
else
    if alpha < thres2
        alpha = thres2;
    end
end

dv = alpha * MinkowskiPooling(D, beta);
end

function D2 = MinkowskiPooling(D, beta)
D2 = mean(D.^beta)^(1/beta);
end