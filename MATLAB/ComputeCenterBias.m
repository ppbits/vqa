% compute center bias map for visual attention
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
