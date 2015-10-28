%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the proposed VQA algorithm
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 17, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Test(database)
if strcmp(database, 'LiveMobile')
    modes = {'phone', 'tablet'};
    for m = 1:length(modes)
        WriteResultLine(modes{m}, 'a');
        for scale = [64 128 256 720]
            for dist_type = [1 3 4 5 -2]
                TestLiveMobile(modes{m}, scale, dist_type);
            end
        end
    end
elseif strcmp(database, 'Live')
    for scale = [64 128 256 432]
        for dist_type = [1 2 3 4 0]
            TestLive(scale, dist_type);
        end
    end
else
    error('Unknown database!')
end
end
