%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the saliency weighted and motion tuned distortion metric
% between two volumes of spatiotemporal energy distributions.
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qPF = MeasureDistortion(volSE1_both, volSE2_both, volSI, applyMotionTuning, applySaliency)
volSE1_m = volSE1_both(:,1); % summed/marginalized motion enery
volSE1 = volSE1_both(:,2); % 4 channel motion enery
clear volSE1_both

volSE2_m = volSE2_both(:,1); % summed/marginalized motion enery
volSE2 = volSE2_both(:,2); % 4 channel motion enery
clear volSE2_both

nChannel = size(volSE1, 1);% number of filters + 1
nFrame = size(volSE1{1,1},3);

N = 3;
filtersPerPlane = N + 1;
nFrequencyPlanes = nChannel-1;
nFilters = filtersPerPlane * nFrequencyPlanes;
tempVolSE = cell(nFilters, 1);
weights4MotionTuning = cell(nFilters,1);
tempVolSE_m = cell(nFrequencyPlanes, 1);
i_filter = 1;
for i = 1:nFrequencyPlanes
    for j = 1:filtersPerPlane
        % difference between rectified filter responses (unmarginalized but normalized) 
        tempVolSE{i_filter} = abs(volSE1{i}{j} - volSE2{i}{j}); 

        % Weight for motion tuning
        weights4MotionTuning{i_filter} = volSE1_m{i};

        i_filter = i_filter + 1;
    end
end

qPF = Minkowski(tempVolSE, 2, volSI, weights4MotionTuning, applyMotionTuning, applySaliency); % motion tuned

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the saliency weighted and motion tuned Minkowski distance
% between two volumes of spatiotemporal energy distributions.
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qPF = Minkowski(absDiffVols, beta, volSI, motionWeights, applyMotionTuning, applySaliency)
nChannel = size(absDiffVols, 1);% number of filters + 1
nPixelPF = size(absDiffVols{1, 1},1) * size(absDiffVols{1, 1}, 2);

if ~applyMotionTuning % No motion tuning
    s4 = absDiffVols{1}.^beta;
    for i = 2:nChannel
        s4 = s4 + absDiffVols{i}.^beta;
    end
else % With motion tuning
    s4 = absDiffVols{1}.^beta.*motionWeights{1};
    for i = 2:nChannel
        s4 = s4 + absDiffVols{i}.^beta.*motionWeights{i};
    end
end

s4 = s4.^(1/beta);
if ~applySaliency
    qPF = sum(sum(s4))/nPixelPF;
else
    s4_SI = s4 .* volSI;% self-info weighted
    volSI_sumPF = sum(sum(volSI));
    qPF = sum(sum(s4_SI))./volSI_sumPF;
end
end
