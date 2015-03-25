% This code tests the algorithm on the Live Mobile VQA Database.
%
% Contact: dante.peng@gmail.com
% Copyright (c) Peng Peng
function TestLiveMobile(dist_type)
if nargin < 1
    dist_type = 0;
end

%% Prepare Data
data_path = 'D:\pengp\LiveMobileVQA';
yuv_path = fullfile(data_path, 'YUV');
mat_path = fullfile(data_path, 'MAT');
frame_size = [720 1280];
% Load 'dist_names', 'org_names', and 'refnames_all'.
load(fullfile(data_path, 'names.mat'));
% Load 'dmos_mobile', 'dmos_tablet', 'std_dmos_mobile', 'std_dmos_tablet'.
load(fullfile(data_path, 'dmos_final.mat'));
% Get distortion types
dist_types = GetDistortionTypes('LiveMobile');
% Create output folder
output_path = 'output';
if(~isdir(output_path))
    fprintf('Creating output folder %s...\n', output_path);
    mkdir(output_path);
end

if dist_type == 0 
	selected_videos = true(length(dmos_mobile),1);
else
    selected_videos = (dist_types == dist_type);
end

ref_videos = refnames_all(selected_videos);
dist_videos = dist_names(selected_videos);
nfile = length(dist_videos);
fprintf('Selected videos: %d.\nDistortion type: %d.\n', nfile, dist_type);

%% Process Data
for scale = [64 128 256]
    for i = 1:nfile
        ref_filename = strcat(ref_videos{i}, '_org');
        fprintf('Reading video %s at scale %d...\n', ref_filename, scale);
        ReadVideo(yuv_path, mat_path, ref_filename, frame_size, scale);
        dist_filename = dist_videos{i};
        fprintf('Reading video %s at scale %d...\n', dist_filename, scale);
        ReadVideo(yuv_path, mat_path, dist_filename, frame_size, scale);
    end
end

fprintf('Done!');
