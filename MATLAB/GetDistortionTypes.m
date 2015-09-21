%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gets a vector representing the distortion types of all distorted videos in a database
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist_types = GetDistortionTypes(database, n)

fprintf('Getting distortion types for the %s database...\n', database);

if(strcmp(database, 'LiveMobile'))
    % Create distion types array
    dist_types = zeros(n, 1);
    res = mod(0:n-1, 20);
    dist_types(res < 4) = 1; % Compression
    dist_types(res >= 4 & res < 8) = 2; % Frame freeze
    dist_types(res >= 8 & res < 11) = 3; % Rate adaptation
    dist_types(res >= 11 & res < 16) = 4; % Temporal dynamics
    dist_types(res >= 16) = 5; % Wireless channel packet-loss
else
    error('Unknown database!');
end

% isCompression = false(200, 1);
% isFreezeStored = false(200, 1);
% isFreezeLive = false(200, 1);
% isRateAdaptaion = false(200, 1);
% isTemporalDynamics = false(200, 1);
% isPacketLoss = false(200, 1);
%
% for i = 0:9
%     isCompression(i*20 + (1:4)) = true;
%     isFreezeStored(i*20 + (5:7)) = true;
%     isFreezeLive(i*20 + 8) = true;
%     isRateAdaptaion(i*20 + (9:11)) = true;
%     isTemporalDynamics(i*20 + (12:16)) = true;
%     isPacketLoss(i*20 + (17:20)) = true;
% end
