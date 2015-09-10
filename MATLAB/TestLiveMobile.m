% This code tests the algorithm on the Live Mobile VQA Database.
%
% Contact: dante.peng@gmail.com
% Copyright (c) Peng Peng
function TestLiveMobile(dist_type, mode)
if nargin < 2
    mode = 'mobile';
end
if nargin < 1
    dist_type = 0;
end

%% Prepare Data
data_path = '/cs/vml2/pengp/LiveMobileVQA';
yuv_path = fullfile(data_path, 'YUV');
mat_path = fullfile(data_path, 'MAT');
score_path = fullfile(data_path, 'scores');

frame_size = [720 1280];
% Load 'dist_names', 'org_names', and 'refnames_all'.
load(fullfile(data_path, 'names.mat'));
% Load 'dmos_mobile', 'dmos_tablet', 'std_dmos_mobile', 'std_dmos_tablet'.
load(fullfile(data_path, 'dmos_final.mat'));
% Get distortion types
dist_types = GetDistortionTypes('LiveMobile');

if dist_type == 0
    selected_videos = true(length(dmos_mobile),1);
else
    selected_videos = (dist_types == dist_type);
end

if strcmp(mode, 'mobile')
    dmos = dmos_mobile(selected_videos);
elseif strcmp(mode, 'tablet')
    dmos = dmos_tablet(selected_videos);
else
    error('Wrong mode. Please select mobile or tablet.');
end
    
ref_filenames = refnames_all(selected_videos);
dist_filenames = dist_names(selected_videos);
nfile = length(dist_filenames);
fprintf('Selected videos: %d.\nDistortion type: %d.\n', nfile, dist_type);

%% Load spatial quality
load(fullfile(score_path, 'liveM_MSSIM.mat'));
spatial_scorePF_all = scoresPerFrameAll(selected_videos);
spatial_mScore_all = 1 - mScoreAll(selected_videos);

%% Compute Overall Quality
spatial_tpScore_all = zeros(nfile, 1);
motion_mScore_all = zeros(nfile, 1);
motion_tpScore_all = zeros(nfile, 1);
column = 2;

for scale = [64 128 256]
    motion_score_folder = fullfile(score_path, int2str(scale));
    if(~isdir(motion_score_folder))
        fprintf('Creating score folder %s...\n', motion_score_folder);
        mkdir(motion_score_folder);
    end
    for i = 1:nfile
        ref_filename = strcat(ref_filenames{i}, '_org');
        dist_filename = dist_filenames{i};
        [mScore, scorePF] = GetMotionScores(motion_score_folder, ref_filename, dist_filename);
        
        motion_mScore_all(i) = mScore(column);
        
        % Temproal pooling
        lambda1 = 1.5; beta = 1; lambda2 = 0.5;
        motion_tpScore_all(i) = temporalPooling(scorePF(:, column), lambda1, beta, lambda2);
        spatial_tpScore_all(i) = temporalPooling(1-spatial_scorePF_all{i}, lambda1, beta, lambda2);
    end
    final_mScore_all = spatial_mScore_all .* motion_mScore_all;
    final_tpScore_all = spatial_tpScore_all .* motion_tpScore_all;
    
    fprintf('Scale = %d\n', scale);
    
    performance(spatial_mScore_all, dmos, true(nfile, 1), 'M->Spatial');
    performance(motion_mScore_all, dmos, true(nfile, 1), 'M->Motion');
    performance(final_mScore_all, dmos, true(nfile, 1), 'M->Product');
    
    performance(spatial_tpScore_all, dmos, true(nfile, 1), 'TP->Spatial');
    performance(motion_tpScore_all, dmos, true(nfile, 1), 'TP->Motion');
    performance(final_tpScore_all, dmos, true(nfile, 1), 'TP->Product');
end

fprintf('\nDone!\n');
