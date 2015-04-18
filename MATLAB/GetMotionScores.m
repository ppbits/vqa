function [mScore, scorePF] = GetMotionScores(score_folder, ref_filename, dist_filename)
score_file = sprintf('%s/%s.mat', score_folder, dist_filename);
if ~exist(score_file, 'file')
    t0 = tic;
    
    fprintf('Reading video %s at scale %d...\n', ref_filename, scale);
    ref_video = ReadVideo(yuv_path, mat_path, ref_filename, frame_size, scale);
    
    fprintf('Reading video %s at scale %d...\n', dist_filename, scale);
    dist_video =ReadVideo(yuv_path, mat_path, dist_filename, frame_size, scale);
    [mScore, scorePF] = AG_MT_SQ_Fast(ref_video, dist_video);
    save(score_file,  'mScore', 'scorePF');
    t1 = toc(t0);
    fprintf('Elapsed time on processing the %d-th video: %3.4f seconds.\n\n', i, t1);
else
    scores_struct = load(score_file);
    mScore = scores_struct.mScore;
    scorePF = scores_struct.scorePF;
end
end