%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute center bias map for visual attention
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 11, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = ComputeCenterBias(h, w)
center = [h/2 w/2];
map  = zeros(h, w);
% normalization factor
D = sqrt(sum(center.^2));
for i = 1:h
    for j = 1:w
        map(i, j) = sqrt(sum(([i j] - center).^2));
    end
end
map = 1-map/D;
end
