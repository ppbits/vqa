%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get quality scores
% 
% Reads the scores from disk if they are already computed and stored. Otherwise, compute the
% quality scores by calling the ComputeQualtyHandle.
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 18, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mScore, scorePF] = GetQualityScores(ComputeQualityHandle, score_folder, yuv_path, mat_path, frame_size, ref_filename, dist_filename, scale)
if(nargin < 8)
    % Full scale (No sub-sampling)
    scale = frame_size(1);
end
score_file = sprintf('%s/%s.mat', score_folder, dist_filename);
if ~exist(score_file, 'file')
    t0 = tic;
    
    fprintf('Reading reference video %s at scale %d...\n', ref_filename, scale);
    ref_video = ReadVideo(yuv_path, mat_path, ref_filename, frame_size, scale);
    
    fprintf('Reading video %s at scale %d...\n', dist_filename, scale);
    dist_video =ReadVideo(yuv_path, mat_path, dist_filename, frame_size, scale);
    [mScore, scorePF] = ComputeQualityHandle(ref_video, dist_video);
    save(score_file,  'mScore', 'scorePF');
    t1 = toc(t0);
    fprintf('Elapsed time on evaluating quality of video %s: %3.4f seconds.\n\n', dist_filename, t1);
else
    fprintf('Reading quality score for video %s at scale %d from file %s\n', dist_filename, scale, score_file);
    scores_struct = load(score_file);
    mScore = scores_struct.mScore;
    scorePF = scores_struct.scorePF;
end
end
