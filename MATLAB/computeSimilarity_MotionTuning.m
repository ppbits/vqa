% This function computes the saliency weighted similarity 
% between two volumes of spatiotemporal energy distributions
%
%% NO Motion Tuning
function [q, qPF] = computeSimilarity_MotionTuning(volSE1_both, volSE2_both, volSI)
volSE1_m = volSE1_both(:,1); % summed/marginalized motion enery
volSE1 = volSE1_both(:,2); % 4 channel motion enery
clear volSE1_both

volSE2_m = volSE2_both(:,1); % summed/marginalized motion enery
volSE2 = volSE2_both(:,2); % 4 channel motion enery
clear volSE2_both

nChannel = size(volSE1, 1);% number of filters + 1
nFrame = size(volSE1{1,1},3);

% motion tuned energy (use only motion information from the reference video)
% for i = 1:nChannel
%     
%     volSE1{i} = volSE1{i} .*  ( volSE1_m{i} );
%     volSE2{i} = volSE2{i} .*  ( volSE1_m{i} );
% 
% % 	volSE1{i} = volSE1{i} .*  (volSE1_m{i} + volSE2_m{i});
% % 	volSE2{i} = volSE2{i} .*  (volSE1_m{i} + volSE2_m{i});
% end
% clear volSE2_m



%% 'Minkowski'
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
        weights4MotionTuning{i_filter} = volSE1_m{i};
        i_filter = i_filter + 1;
    end
    % difference between pure motion descriptor (with marginalization and divisive normalization)
    tempVolSE_m{i} = abs(volSE1_m{i} - volSE2_m{i}); 
end

allBeta = [2];
nBeta = 1;
q_Mink0 = zeros(1,nBeta);
qPF_Mink0 = zeros(nFrame,nBeta);
q_Mink_SI0 = zeros(1,nBeta);
qPF_Mink_SI0 = zeros(nFrame,nBeta);
q_Mink = zeros(1,nBeta);
qPF_Mink = zeros(nFrame,nBeta);
q_Mink_SI = zeros(1,nBeta);
qPF_Mink_SI = zeros(nFrame,nBeta);
q_Mink2 = zeros(1,nBeta);
qPF_Mink2 = zeros(nFrame,nBeta);
q_Mink_SI2 = zeros(1,nBeta);
qPF_Mink_SI2 = zeros(nFrame,nBeta);

for i = 1:nBeta
    beta = allBeta(i);
    [q_Mink0(i), q_Mink_SI0(i), qPF_Mink0(:,i), qPF_Mink_SI0(:,i)] = ...
        Minkowski(tempVolSE, beta, volSI); % dynamics confound with spatial appearance (not tuned to motion)
    [q_Mink(i), q_Mink_SI(i), qPF_Mink(:,i), qPF_Mink_SI(:,i)] = ...
        Minkowski(tempVolSE, beta, volSI, weights4MotionTuning); % motion tuned
    [q_Mink2(i), q_Mink_SI2(i), qPF_Mink2(:,i), qPF_Mink_SI2(:,i)] = ...
        Minkowski(tempVolSE_m, beta, volSI); % pure dynamics
end

% form the output vectors
q = [q_Mink0 q_Mink_SI0 q_Mink q_Mink_SI q_Mink2 q_Mink_SI2];
qPF = [qPF_Mink0 qPF_Mink_SI0 qPF_Mink qPF_Mink_SI qPF_Mink2 qPF_Mink_SI2];
end

function [q, q_SI, qPF, qPF_SI] = Minkowski(absDiffVols, beta, volSI, motionWeights)
volSI_sumPF = sum(sum(volSI));
nChannel = size(absDiffVols, 1);% number of filters + 1
nPixelPF = size(absDiffVols{1, 1},1) * size(absDiffVols{1, 1}, 2);
if nargin < 4
    s4 = absDiffVols{1}.^beta;
    for i = 2:nChannel
        s4 = s4 + absDiffVols{i}.^beta;
    end
    s4 = s4.^(1/beta);
    qPF = sum(sum(s4))/nPixelPF;
    q = mean(qPF);
    s4_SI = s4 .* volSI;% self-info weighted
    qPF_SI = sum(sum(s4_SI))./volSI_sumPF;
    q_SI = mean(qPF_SI);
else
    s4 = absDiffVols{1}.^beta.*motionWeights{1};
    for i = 2:nChannel
        s4 = s4 + absDiffVols{i}.^beta.*motionWeights{i};
    end
    s4 = s4.^(1/beta);
    qPF = sum(sum(s4))/nPixelPF;
    q = mean(qPF);
    s4_SI = s4 .* volSI;% self-info weighted
    qPF_SI = sum(sum(s4_SI))./volSI_sumPF;
    q_SI = mean(qPF_SI);
end

end