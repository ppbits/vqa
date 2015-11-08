%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate on the LIVE Mobile VQA Database
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 11, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TestLiveMobile(mode, scale, dist_type)
if nargin < 3
    dist_type = -2;
end
if nargin < 2
    scale = 64;
end
if nargin < 1
    mode = 'phone';
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

if strcmp(mode, 'phone')
    dmos_all = dmos_mobile;
    dist_names_all = dist_names;
elseif strcmp(mode, 'tablet')
    dmos_all = dmos_tablet;
    dist_names_all = names_tablet;
else
    error('Wrong mode. Please select phone or tablet.');
end

total_file_num = length(dist_names_all);

% Get distortion types
dist_types = GetDistortionTypes('LiveMobile');
 
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
fprintf('Subjective-evaluation device type: %s\n', mode);
fprintf('Number of selected videos: %d\n', nfile);

%% Load spatial quality (Multi-scale SSIM)
%load(fullfile(score_path, 'liveM_MSSIM.mat')); 
%spatial_scorePF_all = scoresPerFrameAll(selected_videos);
%spatial_mScore_all = 1 - mScoreAll(selected_videos);

%% Compute Overall Quality
motion_mScore_all = zeros(nfile, 1);
spatial_mScore_all = zeros(nfile, 1);
motion_tpScore_all = zeros(nfile, 1);
spatial_tpScore_all = zeros(nfile, 1);

motion_score_folder = fullfile(score_path, strcat('motion_', int2str(scale)));
spatial_score_folder = fullfile(score_path, 'spatial');
if(~isdir(motion_score_folder))
    fprintf('Creating score folder %s...\n', motion_score_folder);
    mkdir(motion_score_folder);
end
if(~isdir(spatial_score_folder))
    fprintf('Creating score folder %s...\n', spatial_score_folder);
    mkdir(spatial_score_folder);
end

for i = 1:nfile
    fprintf('Video #: %d\n', i);
    % Get the distorted file name
    dist_filename = dist_filenames{i};
    % Get the original/reference file name based on the distorted file name
    ref_filename = strcat(dist_filenames{i}(1:2), '_org');

    % Get per-frame moiton-quality scores ('motion_scorePF') and the overall motion quality score ('motion_mScore')
    % 'motion_mScore' is the mean value of 'motion_scorePF' over each column
    [motion_mScore, motion_scorePF] = GetQualityScores(@ComputeMotionQuality, motion_score_folder, yuv_path, mat_path, frame_size, ref_filename, dist_filename, scale);
    motion_mScore_all(i) = motion_mScore;
    motion_tpScore_all(i) = TemporalPooling(motion_scorePF);
    
    % Get per-frame spatial-quality scores ('spatial_scorePF') and the overall spatial quality score ('spatial_mScore')
    % 'spatial_mScore' is the mean value of 'spatial_scorePF' over each column
    [spatial_mScore, spatial_scorePF] = GetQualityScores(@ComputeSpatialQuality, spatial_score_folder, yuv_path, mat_path, frame_size, ref_filename, dist_filename);
    spatial_mScore_all(i) = 1 - spatial_mScore;
    spatial_tpScore_all(i) = TemporalPooling(1-spatial_scorePF);
end

% Combine the spatial quality score and motion quality score
final_mScore_all = spatial_mScore_all .* motion_mScore_all;
final_tpScore_all = spatial_tpScore_all .* motion_tpScore_all;
dist_type_name = GetDistortionTypeName('LiveMobile', dist_type);
tag = sprintf('%s - %s\n', mode, dist_type_name);
fprintf(tag);
WriteResult(tag);
Performance(spatial_mScore_all, dmos, true(nfile, 1), strcat('Spatial-', int2str(scale)));
Performance(spatial_tpScore_all, dmos, true(nfile, 1), strcat('TP-Spatial-', int2str(scale)));
Performance(motion_mScore_all, dmos, true(nfile, 1), strcat('Motion-', int2str(scale)));
Performance(motion_tpScore_all, dmos, true(nfile, 1), strcat('TP-Motion-', int2str(scale)));
Performance(final_mScore_all, dmos, true(nfile, 1), strcat('Overall-', int2str(scale)));
Performance(final_tpScore_all, dmos, true(nfile, 1), strcat('TP-Overall-', int2str(scale)));

fprintf('\nDone!\n');