%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the spatial quality of the video using the MS-SSIM index.
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 19, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q, qPF] =  ComputeSpatialQuality(v1, v2)
% addpath(fullfile(pwd, 'lib/ms_ssim'));
addpath(fullfile(pwd, 'lib/msssim'));
nFrame = size(v1, 3);
qPF = zeros(nFrame, 1);

for f = 1 : nFrame
    qPF(f) = msssim(v1(:, :, f), v2(:, :, f));
end
q  = mean(qPF);
end
