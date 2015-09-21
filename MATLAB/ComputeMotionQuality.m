%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attention Guided and Motion Tuned Quality
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q, qPF] =  ComputeMotionQuality(v1, v2, applyMotionTuning, applySaliency)
%% parameter setting
orientationSet = 2;
[~, allOrientations, ~] = SteerableFeatureSetProperties(orientationSet);
G2H2OrG3 = 1; % G3 filter
filter_half_width = 6;
doMarginalize = 2;

nBin = 1000;
nBin_n = 10;
appGaussianFilter = true;
    
nFrame = size(v1,3);
overlap = 0; %default: 2*filter_half_width;
chunkSize =  13;
step = chunkSize - overlap;
nChunk = ceil((nFrame-2*filter_half_width) / step);
if(nChunk <= 0)
    error('nChunk = %d\n', nChunk);
else
    fprintf('nChunk = %d\n', nChunk);
end
nValidScores = nChunk;
qPF = zeros(nValidScores, 1);
scorePos = [0 0];
for c = 1 : nChunk
%    fprintf('Processing Chunk %d\n', c);
    chunkStartPos = 1+(c-1)*step;
    chunkEndPos = chunkStartPos + chunkSize-1;
    if chunkEndPos > nFrame % the last chunk bounded by the last frame of the video 
        chunkEndPos = nFrame;
    end
    chunkRange = chunkStartPos:chunkEndPos;

    % Spatiotemporal Orientation Analysis on the reference video
    [temp1_r, temp2_r] = SpatiotemporalOrientationAnalysis(...
        v1(:,:,chunkRange), allOrientations, ...
        G2H2OrG3, doMarginalize, filter_half_width);

    % Spationtemporal Orentation Analysis on the distorted video
    [~, temp2_d] = SpatiotemporalOrientationAnalysis(...
        v2(:,:,chunkRange), allOrientations, ...
        G2H2OrG3, doMarginalize, filter_half_width);
    
    % Get saliency based on motion descriptor invariant to local luminance contrast
    temp_volSI1 = GetSaliency(temp1_r(:, 1), false, nBin, appGaussianFilter);
    % Get saliency based on motion descriptor confounded by local luminance contrast
    temp_volSI2 = GetSaliency(temp2_r(:, 1), false, nBin_n, appGaussianFilter);
    % Combine the two saliency maps
    temp_volSI4 = CombineSaliency(temp_volSI1, temp_volSI2);
    % Apply center bias
    finalSaliencyWeights = ApplyCenterBias(temp_volSI4);
    clear temp_volSI1  temp_volSI2 temp1_r 
    temp_qPF = MeasureDistortion(temp2_r, temp2_d, finalSaliencyWeights, applyMotionTuning, applySaliency);  
    
    scorePos(1) = scorePos(2) + 1;
    scorePos(2) = scorePos(1) + size(temp_qPF,1)-1;
    scoreRange = scorePos(1):scorePos(2);
    qPF(scoreRange, :) = temp_qPF;
    clear temp2_r temp2_d temp_volSI4
end
q  = mean(qPF);
end
