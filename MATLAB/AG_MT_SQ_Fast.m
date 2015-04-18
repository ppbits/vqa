%% This is the file that contains all the code necessary for
% spatiotemporal orientation analysis
%
% Name: Peng Peng
% Contact: pengp@sfu.ca
% Date: Nov 17, 2012

%% Notes
%% Attention Guided and Motion Tuned Quality

function [q, qPF] =  AG_MT_SQ_Fast(v1, v2, orientationSet)
%% parameter setting
if nargin < 3
   orientationSet = 2;
end

[~, allOrientations, ~] = steerableFeatureSetProperties(orientationSet);
G2H2OrG3 = 1; % G3 filter
filter_half_width = 6;
doMarginalize = 2;

nBin = 1000;
nBin_n = 10;
appGaussianFilter = true;
appCenterBias = false;
    
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
qPF = zeros(nValidScores, 100);
scorePos = [0 0];
for c = 1 : nChunk
%     fprintf('Processing Chunk %d\n', c);
    chunkStartPos = 1+(c-1)*step;
    chunkEndPos = chunkStartPos + chunkSize-1;
    if chunkEndPos > nFrame % the last chunk bounded by the last frame of the video 
        chunkEndPos = nFrame;
    end
    chunkRange = chunkStartPos:chunkEndPos;
    [temp1_r, temp2_r] = spacetimeOrientationAnalysis(...
        v1(:,:,chunkRange), allOrientations, ...
        G2H2OrG3, doMarginalize, filter_half_width);
    [~, temp2_d] = spacetimeOrientationAnalysis(...
        v2(:,:,chunkRange), allOrientations, ...
        G2H2OrG3, doMarginalize, filter_half_width);
    
    temp_volSI1 = getSaliencyVol(temp1_r(:, 1), appCenterBias, nBin, appGaussianFilter);
%     [~, temp_qPF1] = computeSimilarity_MotionTuning(temp2_r, temp2_d, applyCenterBias(temp_volSI1));
    temp_volSI2 = getSaliencyVol(temp2_r(:, 1), appCenterBias, nBin_n, appGaussianFilter);
%     [~, temp_qPF2] = computeSimilarity_MotionTuning(temp2_r, temp2_d, applyCenterBias(temp_volSI2));
%     temp_volSI3 = temp_volSI1 .* temp_volSI2;
%     [~, temp_qPF3] = computeSimilarity_MotionTuning(temp2_r, temp2_d, applyCenterBias(temp_volSI3));
    temp_volSI4 = combineSaliency(temp_volSI1, temp_volSI2);
    clear temp_volSI1  temp_volSI2 temp1_r 
    [~, temp_qPF4] = computeSimilarity_MotionTuning_Fast(temp2_r, temp2_d, applyCenterBias(temp_volSI4));  
%     temp_qPF = [temp_qPF1 temp_qPF2 temp_qPF3 temp_qPF4];
    temp_qPF = temp_qPF4;
    
    scorePos(1) = scorePos(2) + 1;
    scorePos(2) = scorePos(1) + size(temp_qPF,1)-1;
    scoreRange = scorePos(1):scorePos(2);
    qPF(scoreRange,1:size(temp_qPF,2)) = temp_qPF;
    clear temp2_r temp2_d temp_volSI4
end
q  = mean(qPF);
end

% This function computes the saliency weighted similarity 
% between two volumes of spatiotemporal energy distributions
%
%% NO Motion Tuning
function [q, qPF] = computeSimilarity_MotionTuning_Fast(volSE1_both, volSE2_both, volSI)
volSE1_m = volSE1_both(:,1); % summed/marginalized motion enery
volSE1 = volSE1_both(:,2); % 4 channel motion enery
clear volSE1_both

volSE2_m = volSE2_both(:,1); % summed/marginalized motion enery
volSE2 = volSE2_both(:,2); % 4 channel motion enery
clear volSE2_both

nChannel = size(volSE1, 1);% number of filters + 1
nFrame = size(volSE1{1,1},3);

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
%     tempVolSE_m{i} = abs(volSE1_m{i} - volSE2_m{i}); 
end

q_Mink0 = zeros(1,4);
qPF_Mink0 = zeros(nFrame,4);
q_Mink_SI0 = zeros(1,4);
qPF_Mink_SI0 = zeros(nFrame,4);
q_Mink = zeros(1,4);
qPF_Mink = zeros(nFrame,4);
q_Mink_SI = zeros(1,4);
qPF_Mink_SI = zeros(nFrame,4);
q_Mink2 = zeros(1,4);
qPF_Mink2 = zeros(nFrame,4);
q_Mink_SI2 = zeros(1,4);
qPF_Mink_SI2 = zeros(nFrame,4);

for i = 2
%     [q_Mink0(i), q_Mink_SI0(i), qPF_Mink0(:,i), qPF_Mink_SI0(:,i)] = ...
%         Minkowski(tempVolSE, i, volSI); % dynamics confound with spatial appearance (not tuned to motion)
    [q_Mink(i), q_Mink_SI(i), qPF_Mink(:,i), qPF_Mink_SI(:,i)] = ...
        Minkowski(tempVolSE, i, volSI, weights4MotionTuning); % motion tuned
%     [q_Mink2(i), q_Mink_SI2(i), qPF_Mink2(:,i), qPF_Mink_SI2(:,i)] = ...
%         Minkowski(tempVolSE_m, i, volSI); % pure dynamics
end

% form the output vectors
q = q_Mink;
qPF = qPF_Mink;
% q = [q_Mink0 q_Mink_SI0 q_Mink q_Mink_SI q_Mink2 q_Mink_SI2];
% qPF = [qPF_Mink0 qPF_Mink_SI0 qPF_Mink qPF_Mink_SI qPF_Mink2 qPF_Mink_SI2];
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


function [volSE, volSE_n] = spacetimeOrientationAnalysis(videoVol, allOrientations, G2H2OrG3, ...
    doMarginalize, filter_half_width)
if nargin < 5
    filter_half_width = 6;
end

if nargin < 4
    doMarginalize = 1;
end

if nargin < 3
    G2H2OrG3 = 1;           % G3
end


%% define necessary processing parameters
epsilon = 1000; %25 can adjust the eps value (from anywhere from very small (e.g., 1) to very larger (e.g., 10,000) depending on noise in video
padding = 50;       % not used now
agg_size = 5;% -1;      % for computing normalizing energy (5 is normally used)
sampling_rate = 0.5;
subsampleRaw = 0;
% filter_half_width = 6; %default: 6;  % half width of filter (spatially and temporally)
%doMarginalize = 1; % default: 1     % option for specifying whether you want to capture just dynamics (1) or appearance and dynamics (0)
% G2H2OrG3 = 1;           % 0 = G2/H2, 1 = G3
% orientationSet = 2; % we used 2 a lot; 6, used in KTH and UCF sports experiment




%% process each video
clear curVideoVol volG3 volSE tmpVol normVol
clear G2a_img G2b_img G2c_img G2d_img G2e_img G2f_img
clear H2a_img H2b_img H2c_img H2d_img H2e_img H2f_img H2g_img H2h_img H2i_img H2j_img

%% spatially subsample, if desired
if subsampleRaw
    kernel = 1/16*[1 4 6 4 1]';
    curVideoVol = convn(convn(videoVol,...
        kernel, 'same'), kernel', 'same');
    curVideoVol = curVideoVol(1:2:end, 1:2:end, :);
else
    curVideoVol = videoVol;
end

clear videoVol
%% extract oriented energies with subsampling. 0 = G2/H2, 1 = G3
if(G2H2OrG3 == 0)
    fprintf('Initializing the 3D G2 filters ...\n');
    tic;
    [G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img] = imgInit3DG2Sing(curVideoVol);
    fprintf('Done!'); toc
    
    fprintf('Initializing the 3D H2 filters ...\n');
    tic
    [H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img] = imgInit3DH2Sing(curVideoVol);
    fprintf('Done!'); toc
    
    tic
    fprintf('Computing the signatures with G2H2 ...\n');
    [volSE, volSE_n] = compute_signatures_with_G2H2_general(G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img, ...
        H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img, ...
        epsilon, allOrientations);
    fprintf('Done!'); toc
elseif(G2H2OrG3 == 1)
%     fprintf('initializing 3D G3 ...\n');
%     tic;
    volG3 = imgInit3DG3_with_options(curVideoVol,sampling_rate, filter_half_width);
%     fprintf('Done!'); toc
    
%     fprintf('computing signatures with G3 ...\n');
%     tic
    [volSE, volSE_n] = compute_signatures_with_G3_general(volG3, epsilon, allOrientations, agg_size, doMarginalize);
%     fprintf('Done!'); toc;
    % volSE{1}   -- static channel
    % volSE{end} -- epsilon channel
    % sum the responses in the neighbourhood and normalized
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G3   Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G3 = imgInit3DG3_with_options(vol, sampling_rate, filter_half_width)

if (nargin < 2)
    sampling_rate = 0.5;
    filter_half_width = 6;
elseif (nargin < 3)
    filter_half_width = 6;
end

t = sampling_rate*[-filter_half_width:filter_half_width];
C = 0.1840;

f1 = -4*C*t.*(-3 + 2*t.^2).*exp(-t.^2);
f2 = t.*exp(-t.^2);
f3 = -4*C*(-1 + 2*t.^2).*exp(-t.^2);
f4 = exp(-t.^2);
f5 = -8*C*t.*exp(-t.^2);

G3.a = aux_conv3d(vol, f1, f4, f4);
G3.b = aux_conv3d(vol, f3, f2, f4);
G3.c = aux_conv3d(vol, f2, f3, f4);
G3.d = aux_conv3d(vol, f4, f1, f4);
G3.e = aux_conv3d(vol, f3, f4, f2);
G3.f = aux_conv3d(vol, f5, f2, f2);
G3.g = aux_conv3d(vol, f4, f3, f2);
G3.h = aux_conv3d(vol, f2, f4, f3);
G3.i = aux_conv3d(vol, f4, f2, f3);
G3.j = aux_conv3d(vol, f4, f4, f1);

end


%% Called by imgInit3DG3_with_options
function out = aux_conv3d(vol, f1, f2, f3)
filter_length = length(f1);

out = convn(vol,        f1, 'same');   % x-direction
out = convn(out,        f2','same');   % y-direction
out = convn(out,reshape(f3,1,1,filter_length),'valid');   % z-direction

% out = imfilter(vol,        f1,'conv','same','symmetric');   % x-direction
% out = imfilter(out,        f2','conv','same','symmetric');   % y-direction
% out = imfilter(out,reshape(f3,1,1,filter_length),'conv','same','symmetric');   % z-direction
end



%% Function: compute G3 signatures
function [E, E_n] = compute_signatures_with_G3_general(G3, epsilon, orientations, agg_size, doMarginalize)

if ~exist('epsilon', 'var')
    epsilon = 500;
end

if ~exist('orientations', 'var')
    %if nargin < 3
    orientations = [0,  0, 1; % static channel is always the first
        1,  0, 1;
        -1,  0, 1;
        0,  1, 1;
        0, -1, 1;
        0,  1, 0];
end
nOrientations = size(orientations,1)+1;

if ~exist('agg_size', 'var') || isempty(agg_size)
    var_size = 13;
end

if doMarginalize <= 1
    nColumn = 1;
elseif doMarginalize == 2
    nColumn = 2;
end

E = cell(nOrientations, nColumn);

% compute base energies
for e = 1:nOrientations-1
    if(doMarginalize == 1)
        E{e, 1} = imgComputeSpacetimeOrientationEnergyG3(G3, orientations(e,1), orientations(e,2), orientations(e,3), agg_size);
    elseif(doMarginalize == 0)
        E{e, 1} = imgComputeSOEsG3_NoMarginalize(G3, orientations(e,:), agg_size);
    elseif doMarginalize == 2
        [E{e, 1}, E{e,2}]= imgComputeSpacetimeOrientationEnergyG3(G3, orientations(e,1), orientations(e,2), orientations(e,3), agg_size);
    end
end

E{nOrientations, 1} = epsilon*ones(size(E{1,1})); % last orientation is always epsilon
% last orientation is always epsilon
temp = cell(4,1);
for i = 1:4
    temp{i,1} = epsilon*ones(size(E{1,1}));
end
E{nOrientations, 2} = temp;

% normalize, if nt needed
E_n = cell(nOrientations, nColumn);
Esum = E{1, 1};
for e = 2:nOrientations
    Esum = Esum + E{e, 1};
end

for e = 1:nOrientations
    E_n{e, 1} = E{e, 1}./Esum;
end
for e = 1:nOrientations
    for sub_channel = 1:4 % four equally spaced directions on the plane
        E_n{e, 2}{sub_channel,1} = E{e,2}{sub_channel,1} ./ Esum;
    end
end
end


function [motion_energy] = imgComputeSOEsG3_NoMarginalize(G3, orientations, agg_size)

if (nargin < 3 || isempty(agg_size))
    agg_size = 13;
end

g = imgBinomialFilter(agg_size);

% if ~isempty(dir('imgSteer3DG3_fast.mex*'))
%     % use fast addition version
%     motion_energy = ...
%         imgSteer3DG3_fast(single(cos(0)*ua'     + sin(0)*ub'     ), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(1*pi/4)*ua'+ sin(1*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(2*pi/4)*ua'+ sin(2*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(3*pi/4)*ua'+ sin(3*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
% else
motion_energy = imgSteer3DG3(orientations, G3).^2;

%end

% perform aggregation only if kernel is specified
if agg_size > 0
    motion_energy = imfilter(imfilter(imfilter(motion_energy,g),g'),reshape(g,1,1,length(g)));
end

end


function [motion_energy, energy4channel] = imgComputeSpacetimeOrientationEnergyG3(G3, u, v, w, agg_size)
% N+1 equally spaced directions
energy4channel = cell(4,1);

if (nargin < 5 || isempty(agg_size))
    agg_size = 13;
end

g = imgBinomialFilter(agg_size);

n = [u; v; w];
n = n/norm(n);

e1 = [1; 0; 0];
e2 = [0; 1; 0];

if ( abs(acos(dot(n, e1)/norm(n))) > abs(acos(dot(n, e2)/norm(n))) )
    ua = cross(n,e1);
else
    ua = cross(n,e2);
end

ua = ua/norm(ua);

ub = cross(n,ua);
ub = ub/norm(ub);

% if you need more speed, should be able to uncomment this and it will use
% the mexx function (speeds things up roughly by a factor of 2)
if ~isempty(dir('../doesnotexist/imgSteer3DG3_fast.mex*')) % some time get Inf and NAN results
%     fprintf('fast');
    %    use fast addition version    
    energy4channel{1} = imgSteer3DG3_fast(single(cos(0)*ua'     + sin(0)*ub'     ), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    energy4channel{2} = imgSteer3DG3_fast(single(cos(1*pi/4)*ua'+ sin(1*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    energy4channel{3} = imgSteer3DG3_fast(single(cos(2*pi/4)*ua'+ sin(2*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    energy4channel{4} = imgSteer3DG3_fast(single(cos(3*pi/4)*ua'+ sin(3*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    
%     motion_energy = ...
%         imgSteer3DG3_fast(single(cos(0)*ua'     + sin(0)*ub'     ), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(1*pi/4)*ua'+ sin(1*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(2*pi/4)*ua'+ sin(2*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(3*pi/4)*ua'+ sin(3*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
else
%     fprintf('no fast');
    energy4channel{1} =  imgSteer3DG3(cos(0)*ua'      + sin(0)*ub'     , G3).^2 ;
    energy4channel{2} =   imgSteer3DG3(cos(1*pi/4)*ua' + sin(1*pi/4)*ub', G3).^2;
    energy4channel{3} =  imgSteer3DG3(cos(2*pi/4)*ua' + sin(2*pi/4)*ub', G3).^2;
    energy4channel{4} =   imgSteer3DG3(cos(3*pi/4)*ua' + sin(3*pi/4)*ub', G3).^2;
    
%         motion_energy = imgSteer3DG3(cos(0)*ua'      + sin(0)*ub'     , G3).^2 + ...
%         imgSteer3DG3(cos(1*pi/4)*ua' + sin(1*pi/4)*ub', G3).^2 + ...
%         imgSteer3DG3(cos(2*pi/4)*ua' + sin(2*pi/4)*ub', G3).^2 + ...
%         imgSteer3DG3(cos(3*pi/4)*ua' + sin(3*pi/4)*ub', G3).^2;
end

motion_energy = energy4channel{1} + energy4channel{2} + energy4channel{3} + energy4channel{4};

% perform aggregation only if kernel is specified
if agg_size > 0
    motion_energy = imfilter(imfilter(imfilter(motion_energy,g),g'),reshape(g,1,1,length(g)));
    for i = 1:4 
        energy4channel{1}= imfilter(imfilter(imfilter(energy4channel{1},g),g'),reshape(g,1,1,length(g)));
    end
end

end


function [vol] = imgSteer3DG3(direction, G3)

a = direction(1);
b = direction(2);
c = direction(3);

vol =     (a^3)*G3.a ...
    + 3*(a^2)*b*G3.b ...
    + 3*a*(b^2)*G3.c ...
    +     (b^3)*G3.d ...
    + 3*(a^2)*c*G3.e ...
    +   6*a*b*c*G3.f ...
    + 3*(b^2)*c*G3.g ...
    + 3*a*(c^2)*G3.h ...
    + 3*b*(c^2)*G3.i ...
    +     (c^3)*G3.j;
end


% Description: Generates a normalized binomial filter where n represents
% the desired number of taps
function bfilter = imgBinomialFilter(n)

n = n - 1;

bfilter = [];
for i = 0:n
    bfilter = [bfilter nchoosek(n,i)];
end

% --- normalization ---
bfilter = bfilter/sum(bfilter);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2H2 Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E, E_n] = compute_signatures_with_G2H2_general(G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img, ...
    H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img, ...
    epsilon, orientations)

if ~exist('epsilon', 'var')
    epsilon = 500;
end

% add one more orientation for epsilon channel
nOrientations = size(orientations,1)+1;

E = cell(nOrientations, 1);

% compute base energies
for e = 1:nOrientations-1
    E{e} = imgSteer3DG2(orientations(e, :), G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img).^2 + ...
        imgSteer3DH2(orientations(e, :), H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img).^2;
end
E{nOrientations} = epsilon*ones(size(E{1})); % last orientation is always epsilon

% sum across channel to form the denominator
Esum = E{1};
for e = 2:nOrientations
    Esum = Esum + E{e};
end

% normalize by the sum of the consort
E_n = cell(nOrientations, 1);
for e = 1:nOrientations
    E_n{e} = E{e}./Esum;
end

end


function [G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img] = imgInit3DG2(img)
% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---

SAMPLING_RATE = 0.67;
N = (2/sqrt(3))*(2/pi)^(3/4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms t;
% f1 = N*(2*t^2 - 1)*exp(-t^2);
% f2 = exp(-t^2);
% f3 = 2*N*t*exp(-t^2);
% f4 = t*exp(-t^2);
%
% i = SAMPLING_RATE*[-4:4];
% filter_size = length(i);
%
% basis = [subs(f1,t,i); subs(f2,t,i); subs(f3,t,i); subs(f4,t,i)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KC modification to get rid of syms and subs (much faster and works on
% Laplace)
i = SAMPLING_RATE*[-4:4];
filter_size = length(i);

t = i;

f1 = N*(2*t.^2 - 1).*exp(-t.^2);
f2 = exp(-t.^2);
f3 = 2*N*t.*exp(-t.^2);
f4 = t.*exp(-t.^2);
basis = [f1; f2; f3; f4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G2a_img = single(convn(img, reshape(basis(1,:),1,filter_size),'same'));   % x-direction
G2a_img = single(convn(G2a_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
G2a_img = single(convn(G2a_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

G2b_img = single(convn(img, reshape(basis(3,:),1,filter_size),'same'));   % x-direction
G2b_img = single(convn(G2b_img,reshape(basis(4,:),filter_size,1),'same'));   % y-direction
G2b_img = single(convn(G2b_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

G2c_img = single(convn(img, reshape(basis(2,:),1,filter_size),'same'));   % x-direction
G2c_img = single(convn(G2c_img,reshape(basis(1,:),filter_size,1),'same'));   % y-direction
G2c_img = single(convn(G2c_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

G2d_img = single(convn(img, reshape(basis(3,:),1,filter_size),'same'));   % x-direction
G2d_img = single(convn(G2d_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
G2d_img = single(convn(G2d_img,reshape(basis(4,:),1,1,filter_size),'same')); % z-direction

G2e_img = single(convn(img, reshape(basis(2,:),1,filter_size),'same'));   % x-direction
G2e_img = single(convn(G2e_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
G2e_img = single(convn(G2e_img,reshape(basis(4,:),1,1,filter_size),'same')); % z-direction

G2f_img = single(convn(img, reshape(basis(2,:),1,filter_size),'same'));   % x-direction
G2f_img = single(convn(G2f_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
G2f_img = single(convn(G2f_img,reshape(basis(1,:),1,1,filter_size),'same')); % z-direction
end




function [H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img] = imgInit3DH2Sing(img)
% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---

SAMPLING_RATE = 0.67;
N = 0.877776;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms t;
%
% f1 = N*(t^3 - 2.254*t)*exp(-t^2);
% f2 = N*(t^2 - 0.751333)*exp(-t^2);
% f3 = exp(-t^2);
% f4 = N*t*exp(-t^2);
% f5 = t*exp(-t^2);
%
% i = SAMPLING_RATE*[-4:4];
% filter_size = length(i);
%
% basis = [subs(f1,t,i); subs(f2,t,i); subs(f3,t,i); subs(f4,t,i); subs(f5,t,i)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KC modification to get rid of syms and subs
i = SAMPLING_RATE*[-4:4];
filter_size = length(i);

t = i;

f1 = N*(t.^3 - 2.254.*t).*exp(-t.^2);
f2 = N*(t.^2 - 0.751333).*exp(-t.^2);
f3 = exp(-t.^2);
f4 = N*t.*exp(-t.^2);
f5 = t.*exp(-t.^2);

basis = [f1; f2; f3; f4; f5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H2a_img = single(convn(img,    reshape(basis(1,:),1,filter_size),'same'));   % x-direction
H2a_img = single(convn(H2a_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2a_img = single(convn(H2a_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2b_img = single(convn(img,    reshape(basis(2,:),1,filter_size),'same'));   % x-direction
H2b_img = single(convn(H2b_img,reshape(basis(5,:),filter_size,1),'same'));   % y-direction
H2b_img = single(convn(H2b_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2c_img = single(convn(img,    reshape(basis(5,:),1,filter_size),'same'));   % x-direction
H2c_img = single(convn(H2c_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
H2c_img = single(convn(H2c_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2d_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2d_img = single(convn(H2d_img,reshape(basis(1,:),filter_size,1),'same'));   % y-direction
H2d_img = single(convn(H2d_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2e_img = single(convn(img,    reshape(basis(2,:),1,filter_size),'same'));   % x-direction
H2e_img = single(convn(H2e_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2e_img = single(convn(H2e_img,reshape(basis(5,:),1,1,filter_size),'same')); % z-direction

H2f_img = single(convn(img,    reshape(basis(4,:),1,filter_size),'same'));   % x-direction
H2f_img = single(convn(H2f_img,reshape(basis(5,:),filter_size,1),'same'));   % y-direction
H2f_img = single(convn(H2f_img,reshape(basis(5,:),1,1,filter_size),'same')); % z-direction

H2g_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2g_img = single(convn(H2g_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
H2g_img = single(convn(H2g_img,reshape(basis(5,:),1,1,filter_size),'same')); % z-direction

H2h_img = single(convn(img,    reshape(basis(5,:),1,filter_size),'same'));   % x-direction
H2h_img = single(convn(H2h_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2h_img = single(convn(H2h_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

H2i_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2i_img = single(convn(H2i_img,reshape(basis(5,:),filter_size,1),'same'));   % y-direction
H2i_img = single(convn(H2i_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

H2j_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2j_img = single(convn(H2j_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2j_img = single(convn(H2j_img,reshape(basis(1,:),1,1,filter_size),'same')); % z-direction
end



function [img] = imgSteer3DG2(direction, G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img)

% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---

a = direction(1);
b = direction(2);
c = direction(3);

img = (a^2)*G2a_img ...
    + 2*a*b*G2b_img ...
    + (b^2)*G2c_img ...
    + 2*a*c*G2d_img ...
    + 2*b*c*G2e_img ...
    + (c^2)*G2f_img;

end



function [img] = imgSteer3DH2(direction, H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img)

% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---
a = direction(1);
b = direction(2);
c = direction(3);

img = (a^3)*H2a_img ...
    + 3*(a^2)*b*H2b_img ...
    + 3*a*(b^2)*H2c_img ...
    + (b^3)*H2d_img ...
    + 3*(a^2)*c*H2e_img ...
    + 6*a*b*c*H2f_img ...
    + 3*(b^2)*c*H2g_img ...
    + 3*a*(c^2)*H2h_img ...
    + 3*b*(c^2)*H2i_img ...
    + (c^3)*H2j_img;
end






