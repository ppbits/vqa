%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get motion quality scores
% 
% Reads the scores from disk if they are already computed and stored. Otherwise, compute the motion
% quality scores by calling the ComputeMotionQualty method.
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mScore, scorePF] = GetMotionScores(score_folder, yuv_path, mat_path, frame_size, ref_filename, dist_filename, scale)
score_file = sprintf('%s/%s.mat', score_folder, dist_filename);
if ~exist(score_file, 'file')
    t0 = tic;
    
    fprintf('Reading reference video %s at scale %d...\n', ref_filename, scale);
    ref_video = ReadVideo(yuv_path, mat_path, ref_filename, frame_size, scale);
    
    fprintf('Reading video %s at scale %d...\n', dist_filename, scale);
    dist_video =ReadVideo(yuv_path, mat_path, dist_filename, frame_size, scale);
    [mScore, scorePF] = ComputeMotionQuality(ref_video, dist_video, true, true); %false);
    save(score_file,  'mScore', 'scorePF');
    t1 = toc(t0);
    fprintf('Elapsed time on processing the %d-th video: %3.4f seconds.\n\n', i, t1);
else
    scores_struct = load(score_file);
    mScore = scores_struct.mScore;
    scorePF = scores_struct.scorePF;
end
end
