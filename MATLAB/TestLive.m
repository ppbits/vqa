%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate on the LIVE VQA Database
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 18, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TestLive(scale, dist_type)
if nargin < 2
    scale = 64;
end
if nargin < 1
    dist_type = -2;
end

%% Prepare Data
if isunix
    data_path = '/cs/vml2/pengp/LiveVQA';
else
    data_path = '\\bluebell\vml2\pengp\LiveVQA';
end
    
yuv_path = fullfile(data_path, 'YUV');
mat_path = fullfile(data_path, 'MAT');
score_path = fullfile(data_path, 'scores');

frame_size = [768 432];
org_names = {'pa', 'rb', 'rh', 'tr', 'st', 'sf', 'bs', 'sh', 'mc', 'pr'};
frame_rates = [25 25 25 25 25 25 25 50 50 50];
% Load the filenames of the distorted videos
dist_names_all = importdata(fullfile(data_path, 'live_video_quality_seqs.txt'));
total_file_num = length(dist_names_all);
subjective_data = importdata(fullfile(data_path, 'live_video_quality_data.txt'));
dmos_all = subjective_data(:, 1);

% Get distortion types
dist_types = GetDistortionTypes('Live');
 
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
% TODO: To be calculated
load(fullfile(score_path, 'live_MSSIM.mat')); 
spatial_scorePF_all = scoresPerFrameAll(selected_videos);
spatial_mScore_all = 1 - mScoreAll(selected_videos);

%% Compute Overall Quality
motion_mScore_all = zeros(nfile, 1);

motion_score_folder = fullfile(score_path, strcat('motion_', int2str(scale)));
if(~isdir(motion_score_folder))
    fprintf('Creating score folder %s...\n', motion_score_folder);
    mkdir(motion_score_folder);
end
for i = 1:nfile
    fprintf('Video #: %d\n', i);
    % Get the distorted file name
    dist_filename = dist_filenames{i};
    % Get the original/reference file name based on the distorted file name
    parts = strsplit(dist_filename, '_');
    ref_filename = strcat(parts{1}(1:2), '1_', parts{2});

    % Get per-frame moiton-quality scores ('scorePF') and the overall quality score ('mScore')
    % 'mScore' is the mean value of 'scorePF' over each column
    [mScore, scorePF] = GetQualityScores(@ComputeMotionQuality, motion_score_folder, yuv_path, mat_path, frame_size, ref_filename, dist_filename, scale);

    % Select the metric
    motion_mScore_all(i) = mScore;
end

% Combine the spatial quality score and motion quality score
final_mScore_all = spatial_mScore_all .* motion_mScore_all;
dist_type_name = GetDistortionTypeName('Live', dist_type);
tag = sprintf('LIVE - %s - %d\n', dist_type_name, scale);
fprintf(tag);
WriteResult(tag);
Performance(spatial_mScore_all, dmos, true(nfile, 1), strcat('Spatial-', int2str(scale)));
Performance(motion_mScore_all, dmos, true(nfile, 1), strcat('Motion-', int2str(scale)));
Performance(final_mScore_all, dmos, true(nfile, 1), strcat('Overall-', int2str(scale)));

fprintf('\nDone!\n');
