function vol = applyCenterBias(volSI)
centerBiasMap = computeCenterBias(size(volSI,1), size(volSI,2));
cbWeights = repmat(centerBiasMap, [1,1,size(volSI,3)]);
vol = volSI .* cbWeights;
end