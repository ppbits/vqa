function volSI =  getSaliencyWeights(volSE, volSE_n)
% fprintf('computing self-information map...\n');	
% t_si_start = tic;

nBin = 1000; 
nBin_n = 10;
appGaussianFilter = true;
appCenterBias = true;
volSI = getSaliencyVol(volSE(:, 2), appCenterBias, nBin, appGaussianFilter);
% volSI_n = getSaliencyVol(volSE_n(:, 2), appCenterBias, nBin_n, appGaussianFilter);
% volSI = combineSaliency(volSI_n, volSI);
% volSI = applyCenterBias(volSI);
% volSI = volSI .* volSI_n;

% volSI = applyCenterBias(volSI .* volSI_n);
% t_si_end = toc(t_si_start);
% fprintf('Done! Elapsed time %3.4f seconds.\n', t_si_end);
end

function s = combineSaliency(s1, s2)
len = size(s1,3);

s = zeros(size(s1));
for i = 1:len
    f1 = s1(:,:,i);
    f2 = s2(:,:,i);
    min1 = min(f1(:));
    min2 = min(f2(:));
    max1 = max(f1(:));
    max2 = max(f2(:));
    f1_n = (f1 - min1)/(max1-min1);
    f2_n = (f2 - min2)/(max2-min2);
    
    
%     std1 = std(f1(:));
%     std2 = std(f2(:));
    s(:,:,i) = f1_n + f1_n.*f2_n;
end

end