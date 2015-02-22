% interface.m -- interface between your M-file and matlab_script
%
% Modified from http://scv.bu.edu/computation/bladecenter/multiple_matlab_tasks.html
%
function interface_liveBigFeaSet(o_i)
if nargin < 1
    o_i = 0;
end
%o_i = 0;

if ~isunix
    addpath(genpath('Z:\IQA/Toolbox/metrix_mux'));
    videoPath = '\\bluebell\vml2\pengp/liveV/vqdatabase/videos/MAT128';
    rootPath = '\\bluebell\vml2\pengp/liveV/vqdatabase';
    siPath = '\\bluebell\vml2\pengp/liveV/vqdatabase/self-info';
else
    addpath(genpath('/cs/grad1/pengp/IQA/Toolbox/metrix_mux'));
    videoPath = '/cs/vml2/pengp/liveV/vqdatabase/videos/MAT128';
    rootPath = '/cs/vml2/pengp/liveV/vqdatabase';
    siPath = '/cs/vml2/pengp/liveV/vqdatabase/self-info';
end

resPath  = sprintf('%s/../Nov26/128_PureMotion_Set3', rootPath);
if ~isdir(resPath)
    mkdir(resPath);
end

load(sprintf('%s/dmos.mat', rootPath));
load(sprintf('%s/dist_types.mat', rootPath));
load(sprintf('%s/ref_names.mat', rootPath));
load(sprintf('%s/dist_names.mat', rootPath));
if o_i == 0
    selected_videos = (dist_types > 0);
else
    selected_videos = (dist_types == o_i);
end
ref_videos = ref_names(selected_videos);
dist_videos = dist_names(selected_videos);
nfile = length(dist_videos);
disp(nfile)


for i = nfile:-1:1
    t0 = tic;
    disp(i)
    fprintf('evaluating the %d-th video ...\n', i);
    ref_filename = ref_videos{i};
    dist_filename = dist_videos{i};
    
    resFile = sprintf('%s/%s.mat', resPath, dist_filename(1:end-4));
    if exist(resFile, 'file') > 0
        fprintf('%s exists!\n', resFile);
        continue;
    else
        fprintf('%s does not exist!\n', resFile);
    end
    clear v1 v2
    load([videoPath '/' ref_filename(1:end-4) '.mat']);
    v1 = videoVol;clear videoVol
    load([videoPath '/' dist_filename(1:end-4) '.mat']);
    v2 = videoVol;clear videoVol
    

    % Propose method
    G2H2OrG3 = 1; % G3 filter
    [mScore, scorePF]    = AttentionGuidedMotionTunedVQA(v1, v2, 3);
    save(resFile,  'mScore', 'scorePF');
    t1 = toc(t0);
    fprintf('Elapsed time on processing the %d-th video: %3.4f seconds.\n', i, t1);
end


