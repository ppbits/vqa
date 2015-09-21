function volSI = GetSaliency(volSE, useCenterBias, nbin, gfilter)
if nargin < 4
    gfilter = false;
end

if nargin < 3
    nbin = 100;
end

if nargin < 2
	useCenterBias = false;
end

volSI = ComputeSelfInformation(volSE, nbin, gfilter);

if useCenterBias
    volSI = ApplyCenterBias(volSI);
end

end
