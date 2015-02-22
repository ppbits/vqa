function volSI = computeSelfInfo(volSE, nbin, gfilter)
if nargin < 3 % Apply Gaussian filter
    gfilter = true;
end
if nargin < 2
    nbin = 100;
end
nFeature = size(volSE, 1);
nFrame  = size(volSE{1}, 3);
h = size(volSE{1}, 1);
w = size(volSE{1}, 2);
volSI = zeros(h, w, nFrame);

weights = ones(100,1);
% weights(5:2:11) = 2;
% weights(12) = 0;

joint_p = ones(h*w, nFrame);
for i = 1: nFeature-1
    for j = 1:nFrame
        frame = volSE{i}(:,:,j);
        frame = frame(:);
        [n, bin] = densityEstimate(frame, nbin);
        p = n/w/h;% probability distribution
        pixel_p = p(bin); % probability of the value at each pixel
        joint_p(:, j) = weights(i) * joint_p(:, j) .* pixel_p; % joint probability over all features (at each pixel)
    end 
end

G = fspecial('gaussian', [10 10], 5);
for j = 1:nFrame
    selfInfo = -log2(joint_p(:, j));
    selfInfoMap = reshape(selfInfo, h, w);
    if gfilter
        volSI(:,:,j) =  imfilter(selfInfoMap, G, 'same');
    else
        volSI(:,:,j) = selfInfoMap;
    end
end


end

% 1000 bin histogram density estimate
function [n, bin] = densityEstimate(img, nbin)
minVal = min(img(:));
maxVal  = max(img(:));
% [f, xi] = ksdensity(img(:));
edges = linspace(minVal, maxVal, nbin);
[n, bin] = histc(img, edges);
end
