function volSI = getSaliencyVol(volSE, useCenterBias, nbin, gfilter)
if nargin < 4
    gfilter = false;
end

if nargin < 3
    nbin = 100;
end

if nargin < 2
	useCenterBias = false;
end

volSI = computeSelfInfo(volSE, nbin, gfilter);

if useCenterBias
	centerBiasMap = computeCenterBias(size(volSI,1), size(volSI,2));
	cbWeights = repmat(centerBiasMap, [1,1,size(volSI,3)]);
	volSI = volSI .* cbWeights;
end

end
