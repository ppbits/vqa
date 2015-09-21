%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply center bias to a saliency map
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vol = ApplyCenterBias(volSI)
centerBiasMap = ComputeCenterBias(size(volSI,1), size(volSI,2));
cbWeights = repmat(centerBiasMap, [1,1,size(volSI,3)]);
vol = volSI .* cbWeights;
end
