% This function computes the saliency weighted similarity 
% between two volumes of spatiotemporal energy distributions
%
%% NO Motion Tuning
function [q, qPF] = computeSimilarity2(volSE1_both, volSE2_both, volSI)
volSE1 = volSE1_both(:,1);
volSE2 = volSE2_both(:,1);
volSE1_m = volSE1_both(:,2);
volSE2_m = volSE2_both(:,2);
clear volSE1_both 

volSI_sumPF = sum(sum(volSI));
nChannel = size(volSE1, 1);% number of filters + 1
nPixelPF = size(volSE1{1, 1},1) * size(volSE1{1, 1}, 2);
nFrame = size(volSE1{1,1},3);
% motion tuned energy (use only motion information from the reference video)
% for i = 1:nChannel
%     
%     volSE1{i} = volSE1{i} .*  ( volSE1_m{i} );
%     volSE2{i} = volSE2{i} .*  ( volSE1_m{i} );
% 
% % 	volSE1{i} = volSE1{i} .*  (volSE1_m{i} + volSE2_m{i});
% % 	volSE2{i} = volSE2{i} .*  (volSE1_m{i} + volSE2_m{i});
% end
% clear volSE2_m

%% 'Bhattacharyya'
tempVolSE = cell(nChannel, 1);
for i = 1:nChannel
    tempVolSE{i} = volSE1{i}.*volSE2{i};
    tempVolSE{i} = sqrt(tempVolSE{i});
end
s = tempVolSE{1}.*volSE1_m{1};
for i =2:nChannel % use the epsilon channel
    s = s + tempVolSE{i}.*volSE1_m{i};
end
qPF_Bhattacharyya = sum(sum(s))/nPixelPF;
q_Bhattacharyya = mean(qPF_Bhattacharyya);
s_SI = s .* volSI;
qPF_Bhattacharyya_SI = sum(sum(s_SI))./volSI_sumPF;
q_Bhattacharyya_SI = mean(qPF_Bhattacharyya_SI);
clear tempVolSE s


tempVolSE = cell(nChannel, 1);
for i = 1:nChannel
    tempVolSE{i} = volSE1_m{i}.*volSE2_m{i};
    tempVolSE{i} = sqrt(tempVolSE{i});
end
s = tempVolSE{1};
for i =2:nChannel % use the epsilon channel
    s = s + tempVolSE{i};
end
qPF_Bhattacharyya2 = sum(sum(s))/nPixelPF;
q_Bhattacharyya2 = mean(qPF_Bhattacharyya2);
s_SI = s .* volSI;
qPF_Bhattacharyya_SI2 = sum(sum(s_SI))./volSI_sumPF;
q_Bhattacharyya_SI2 = mean(qPF_Bhattacharyya_SI2);
qPF_Bhattacharyya = reshape(qPF_Bhattacharyya(1,1,:), nFrame, 1);
qPF_Bhattacharyya2 = reshape(qPF_Bhattacharyya2(1,1,:), nFrame, 1);
qPF_Bhattacharyya_SI = reshape(qPF_Bhattacharyya_SI(1,1,:), nFrame, 1);
qPF_Bhattacharyya_SI2 = reshape(qPF_Bhattacharyya_SI2(1,1,:), nFrame, 1);
clear tempVolSE

%% 'Minkowski'
tempVolSE = cell(nChannel-1, 1);
tempVolSE_m = cell(nChannel-1, 1);

for i = 1:nChannel-1
    tempVolSE{i} = abs(volSE1{i} - volSE2{i});
    tempVolSE_m{i} = abs(volSE1_m{i} - volSE2_m{i});
end


q_Mink = zeros(1,4);
qPF_Mink = zeros(nFrame,4);
q_Mink_SI = zeros(1,4);
qPF_Mink_SI = zeros(nFrame,4);
q_Mink2 = zeros(1,4);
qPF_Mink2 = zeros(nFrame,4);
q_Mink_SI2 = zeros(1,4);
qPF_Mink_SI2 = zeros(nFrame,4);

for i = [1 2 3 4]
    [q_Mink(i), q_Mink_SI(i), qPF_Mink(:,i), qPF_Mink_SI(:,i)] = ...
        Minkowski(tempVolSE, i, volSI, volSE1_m); % motion tuned
    [q_Mink2(i), q_Mink_SI2(i), qPF_Mink2(:,i), qPF_Mink_SI2(:,i)] = ...
        Minkowski(tempVolSE_m, i, volSI); % pure dynamics
end

% form the output vectors
q = [q_Bhattacharyya q_Bhattacharyya_SI q_Bhattacharyya2 q_Bhattacharyya_SI2  q_Mink q_Mink_SI q_Mink2 q_Mink_SI2];
qPF = [qPF_Bhattacharyya qPF_Bhattacharyya_SI qPF_Bhattacharyya2 qPF_Bhattacharyya_SI2 qPF_Mink qPF_Mink_SI qPF_Mink2 qPF_Mink_SI2];
end

function [q, q_SI, qPF, qPF_SI] = Minkowski(absDiffVols, beta, volSI, motionWeights)
volSI_sumPF = sum(sum(volSI));
nChannel = size(absDiffVols, 1);% number of filters + 1
nPixelPF = size(absDiffVols{1, 1},1) * size(absDiffVols{1, 1}, 2);
if nargin < 4
    s4 = absDiffVols{1}.^beta;
    for i = 2:nChannel
        s4 = s4 + absDiffVols{i}.^beta;
    end
    s4 = s4.^(1/beta);
    qPF = sum(sum(s4))/nPixelPF;
    q = mean(qPF);
    s4_SI = s4 .* volSI;% self-info weighted
    qPF_SI = sum(sum(s4_SI))./volSI_sumPF;
    q_SI = mean(qPF_SI);
else
    s4 = absDiffVols{1}.^beta.*motionWeights{1};
    for i = 2:nChannel
        s4 = s4 + absDiffVols{i}.^beta.*motionWeights{i};
    end
    s4 = s4.^(1/beta);
    qPF = sum(sum(s4))/nPixelPF;
    q = mean(qPF);
    s4_SI = s4 .* volSI;% self-info weighted
    qPF_SI = sum(sum(s4_SI))./volSI_sumPF;
    q_SI = mean(qPF_SI);
end

end