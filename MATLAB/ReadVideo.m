%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads a video at a specific scale.
% 
% If the video has already been scaled and stored on the disk as .mat file, read in the .mat file.
% Otherwise, reads in the .yuv file, down-sample it to the desired scale and store it as .mat file.
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function video = ReadVideo(yuv_path, mat_path, filename, size, scale)
    mat_folder = sprintf('%s/%d', mat_path, scale);
    if ~isdir(mat_folder)
        sprintf('Creating output folder %s\n', mat_folder);
        mkdir(mat_folder);
    end
    yuv_file = sprintf('%s/%s.yuv', yuv_path, filename);
    mat_file = sprintf('%s/%s.mat', mat_folder, filename);
    
    if ~exist(mat_file, 'file')
        fprintf('Reading YUV file %s...\n', yuv_file);
        video = ReadYUV(yuv_file, size);
        video = DownSample(video, scale);
        fprintf('Saving video to %s...\n', mat_file);
        save(mat_file, 'video');
    else
        fprintf('Loading video from %s...\n', mat_file);
        video_struct = load(mat_file);
        video = video_struct.video;
    end
end
