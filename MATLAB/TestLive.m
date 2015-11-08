%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate on the LIVE VQA Database
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 18, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TestLive(scale, dist_type)
if nargin < 2
    dist_type = 0;
end
if nargin < 1
    scale = 64;
end

%% Prepare Data
if isunix
    data_path = '/cs/vml2/pengp/LiveVQA_1';
else
    data_path = '\\bluebell\vml2\pengp\LiveVQA_1';
end
    
yuv_path = fullfile(data_path, 'YUV');
mat_path = fullfile(data_path, 'MAT');
score_path = fullfile(data_path, 'scores');

frame_size = [432 768];
%org_names = {'pa', 'rb', 'rh', 'tr', 'st', 'sf', 'bs', 'sh', 'mc', 'pr'};
%frame_rates = [25 25 25 25 25 25 25 50 50 50];

% Load the filenames of the distorted videos
dist_names_all = importdata(fullfile(data_path, 'live_video_quality_seqs.txt'));
total_file_num = length(dist_names_all);
for i=1:total_file_num
dist_names_all{i}=dist_names_all{i}(1:end-4);
end
subjective_data = importdata(fullfile(data_path, 'live_video_quality_data.txt'));
dmos_all = subjective_data(:, 1);

% Get distortion types
dist_types = GetDistortionTypes('Live', total_file_num);
 
if dist_type == 0 % all distortion types
    selected_videos = true(total_file_num,1);
elseif dist_type > 0 % select a single distortion
    selected_videos = (dist_types == dist_type);
end

dmos = dmos_all(selected_videos);
dist_filenames = dist_names_all(selected_videos);
nfile = length(dist_filenames);
%fprintf('Subjective-evaluation device type: %s\n', mode);
fprintf('Number of selected videos: %d\n', nfile);

%% Load spatial quality (Multi-scale SSIM)
% TODO: To be calculated
%load(fullfile(data_path,'spatial_mScore_all.mat'));
%spatial_mScore_all = spatial_mScore_all(selected_videos);

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
    parts = strsplit(dist_filename, '_');
    ref_filename = strcat(parts{1}(1:2), '1_', parts{2});

    % Get per-frame moiton-quality scores ('scorePF') and the overall quality score ('mScore')
    % 'mScore' is the mean value of 'scorePF' over each column
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
dist_type_name = GetDistortionTypeName('Live', dist_type);
tag = sprintf('LIVE - %s - %d\n', dist_type_name, scale);
fprintf(tag);
WriteResult(tag);
Performance(spatial_mScore_all, dmos, true(nfile, 1), strcat('Spatial-', int2str(scale)));
Performance(spatial_tpScore_all, dmos, true(nfile, 1), strcat('TP-Spatial-', int2str(scale)));
Performance(motion_mScore_all, dmos, true(nfile, 1), strcat('Motion-', int2str(scale)));
Performance(motion_tpScore_all, dmos, true(nfile, 1), strcat('TP-Motion-', int2str(scale)));
Performance(final_mScore_all, dmos, true(nfile, 1), strcat('Overall-', int2str(scale)));
Performance(final_tpScore_all, dmos, true(nfile, 1), strcat('TP-Overall-', int2str(scale)));
fprintf('\nDone!\n');
