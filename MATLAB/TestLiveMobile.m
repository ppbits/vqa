%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate on the LIVE Mobile VQA Database
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TestLiveMobile(dist_type, mode, scale)
if nargin < 3
    scale = 64;
end
if nargin < 2
    mode = 'mobile';
end
if nargin < 1
    dist_type = 0;
end

%% Prepare Data
if isunix
    data_path = '/cs/vml2/pengp/LiveMobileVQA';
else
    data_path = '\\bluebell\vml2\pengp\LiveMobileVQA';
end
    
yuv_path = fullfile(data_path, 'YUV');
mat_path = fullfile(data_path, 'MAT');
score_path = fullfile(data_path, 'scores');

frame_size = [720 1280];

% Load 'dist_names', 'org_names', and 'refnames_all', 'names_tablet'.
load(fullfile(data_path, 'names.mat'));
load(fullfile(data_path, 'names_tablet.mat'));

% Load 'dmos_mobile', 'dmos_tablet', 'std_dmos_mobile', 'std_dmos_tablet'.
load(fullfile(data_path, 'dmos_final.mat'));

if strcmp(mode, 'mobile')
    dmos_all = dmos_mobile;
    dist_names_all = dist_names;
elseif strcmp(mode, 'tablet')
    dmos_all = dmos_tablet;
    dist_names_all = names_tablet;
else
    error('Wrong mode. Please select mobile or tablet.');
end

total_file_num = length(dist_names_all);

% Get distortion types
dist_types = GetDistortionTypes('LiveMobile', total_file_num);
 
if dist_type == 0 % all distortion types
    selected_videos = true(total_file_num,1);
elseif dist_type > 0 % select a single distortion
    selected_videos = (dist_types == dist_type);
else % exclude the distortion type 'abs(dist_type)'
    selected_videos = (dist_types ~= abs(dist_type));
end

dmos = dmos_all(selected_videos);
dist_filenames = dist_names(selected_videos);
nfile = length(dist_filenames);
fprintf('Selected videos: %d.\nDistortion type: %d.\n', nfile, dist_type);

%% Load spatial quality (Multi-scale SSIM)
load(fullfile(score_path, 'liveM_MSSIM.mat')); 
spatial_scorePF_all = scoresPerFrameAll(selected_videos);
spatial_mScore_all = 1 - mScoreAll(selected_videos);

%% Compute Overall Quality
motion_mScore_all = zeros(nfile, 1);

motion_score_folder = fullfile(score_path, int2str(scale));
if(~isdir(motion_score_folder))
    fprintf('Creating score folder %s...\n', motion_score_folder);
    mkdir(motion_score_folder);
end
for i = 1:nfile
    % Get the distorted file name
    dist_filename = dist_filenames{i};
    % Get the original/reference file name based on the distorted file name
    ref_filename = strcat(dist_filenames{i}(1:2), '_org');

    % Get per-frame moiton-quality scores ('scorePF') and the overall quality score ('mScore')
    % 'mScore' is the mean value of 'scorePF' over each column
    [mScore, scorePF] = GetMotionScores(motion_score_folder, yuv_path, mat_path, frame_size, ref_filename, dist_filename, scale);

    disp(size(mScore));
    % Select the metric
    motion_mScore_all(i) = mScore;
end

% Combine the spatial quality score and motion quality score
final_mScore_all = spatial_mScore_all .* motion_mScore_all;

fprintf('Scale = %d\n', scale);
Performance(spatial_mScore_all, dmos, true(nfile, 1), strcat('Spatial-', int2str(scale)));
Performance(motion_mScore_all, dmos, true(nfile, 1), strcat('Motion-', int2str(scale)));
Performance(final_mScore_all, dmos, true(nfile, 1), strcat('Overall-', int2str(scale)));

fprintf('\nDone!\n');
