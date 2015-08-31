function [srcc, krcc, plcc, RMSE] = performance(predictions, dmos, selected_videos, tag, printToFile)
if(argin < 5>)
    printToFile = false;
end
if nargin < 4
    tag = 'Performance';
end
if nargin < 3
    selected_videos = true(length(dmos), 1);
end
predictions = predictions(selected_videos);
dmos = dmos(selected_videos);
srcc = corr(predictions, dmos, 'type', 'Spearman');
krcc = corr(predictions, dmos, 'type', 'Kendall');

%  plcc = 0;
%  RMSE = 0;

% %    plcc = corr(predictions(selected_videos), dmos(selected_videos), 'type', 'Pearson');

mos = dmos;
beta(1) = max(mos);
beta(2) = min(mos);
beta(3) = mean(predictions);
beta(4) = 0.1;
beta(5) = 0.1;
%fitting a curve using the data
[bayta ehat,J] = nlinfit(predictions,mos,@logistic,beta, statset('Display','off'));
%given a ssim value, predict the correspoing mos (ypre) using the fitted curve
[ypre junk] = nlpredci(@logistic,predictions,bayta,ehat,J);
RMSE = sqrt(sum((ypre - mos).^2) / length(mos));%root meas squared error
plcc = corr(mos, ypre, 'type','Pearson'); %pearson linear coefficient

% 
fprintf([tag ':\tSRCC = %3.4f, KRCC = %3.4f, PLCC = %3.4f, RMSE = %3.4f\n'], ...
    abs(srcc), abs(krcc), abs(plcc), RMSE);
%

if printToFile
    fid = fopen('Results.txt', 'a');
    if fid~=0
        fprintf(fid, [tag ':\tSRCC = %3.4f, KRCC = %3.4f, PLCC = %3.4f, RMSE = %3.4f\r\n'], ...
            abs(srcc), abs(krcc), abs(plcc), abs(RMSE));
    
    end
    fprintf(fid, '\r\n');
    fclose(fid);
end