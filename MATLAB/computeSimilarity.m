% This function computes the saliency weighted similarity 
% between two volumes of spatiotemporal energy distributions
%
function [q, qPF] = computeSimilarity(volSE1_both, volSE2, volSI)
volSE1 = volSE1_both(:,1);
volSE1_m = volSE1_both(:,2);
% volSE2_m = volSE2_both(:,2);
clear volSE1_both 

volSI_sumPF = sum(sum(volSI));
nChannel = size(volSE1, 1);% number of filters + 1
nPixelPF = size(volSE1{1, 1},1) * size(volSE1{1, 1}, 2);

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

clear tempVolSE

%% 'Minkowski'
tempVolSE = cell(nChannel-1, 1);
for i = 1:nChannel-1
    tempVolSE{i} = volSE1{i} - volSE2{i};
end

% Manhattan distance
s1 = abs(tempVolSE{1}).*volSE1_m{1};
for i = 2:nChannel-1
    s1 = s1 + abs(tempVolSE{i}).*volSE1_m{i}; 
end
qPF_Manhattan = sum(sum(s1))/nPixelPF;
q_Manhattan = mean(qPF_Manhattan);
s1_SI = s1 .* volSI; % self-info weighted
qPF_Manhattan_SI = sum(sum(s1_SI))./volSI_sumPF;
q_Manhattan_SI = mean(qPF_Manhattan_SI);
clear s1

% Euclidean distance
s2 = tempVolSE{1}.^2.*volSE1_m{1};
for i = 2:nChannel-1
    s2 = s2 + tempVolSE{i}.^2.*volSE1_m{i};
end
s2 = sqrt(s2);
qPF_Euclidean = sum(sum(s2))/nPixelPF;
q_Euclidean = mean(qPF_Euclidean);
s2_SI = s2 .* volSI;% self-info weighted
qPF_Euclidean_SI = sum(sum(s2_SI))./volSI_sumPF;
q_Euclidean_SI = mean(qPF_Euclidean_SI);
clear s2 

% Minkowski distance with beta = 4
beta = 4;
s4 = tempVolSE{1}.^beta.*volSE1_m{1};
for i = 2:nChannel-1
    s4 = s4 + tempVolSE{i}.^beta.*volSE1_m{i};
end
s4 = s4.^(1/beta);
qPF_Minkowski4 = sum(sum(s4))/nPixelPF;
q_Minkowski4 = mean(qPF_Minkowski4);
s4_SI = s4 .* volSI;% self-info weighted
qPF_Minkowski4_SI = sum(sum(s4_SI))./volSI_sumPF;
q_Minkowski4_SI = mean(qPF_Minkowski4_SI);
clear s4 

% form the output vectors
qPF = [qPF_Bhattacharyya qPF_Manhattan qPF_Euclidean qPF_Minkowski4 qPF_Bhattacharyya_SI  qPF_Manhattan_SI  qPF_Euclidean_SI qPF_Minkowski4_SI];
q = [q_Bhattacharyya q_Manhattan  q_Euclidean q_Minkowski4 q_Bhattacharyya_SI q_Manhattan_SI  q_Euclidean_SI q_Minkowski4_SI];
end
