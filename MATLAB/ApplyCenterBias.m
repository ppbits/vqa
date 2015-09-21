function vol = ApplyCenterBias(volSI)
centerBiasMap = ComputeCenterBias(size(volSI,1), size(volSI,2));
cbWeights = repmat(centerBiasMap, [1,1,size(volSI,3)]);
vol = volSI .* cbWeights;
end
