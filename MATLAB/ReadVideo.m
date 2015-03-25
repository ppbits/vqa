
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
        video = Read_YUV(yuv_file, size);
        video = DownSample(video, scale);
        fprintf('Saving video to %s...\n', mat_file);
        save(mat_file, 'video');
    else
        fprintf('Loading video from %s...\n', mat_file);
        video = load(mat_file);
    end
end