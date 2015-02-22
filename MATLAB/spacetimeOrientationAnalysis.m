function [volSE, volSE_n] = spacetimeOrientationAnalysis(videoVol, G2H2OrG3, orientationSet, doMarginalize)
addpath(genpath('./EnergyCode'));

if nargin < 4
    doMarginalize = 1;
end

if nargin < 3
	orientationSet = 2;
end

if nargin < 2
    G2H2OrG3 = 0;           % G2H2
end

%% define necessary processing parameters
epsilon = 1000; %25 can adjust the eps value (from anywhere from very small (e.g., 1) to very larger (e.g., 10,000) depending on noise in video
padding = 50;       % not used now
agg_size = 5; %-1;      % for computing normalizing energy (5 is normally used)
sampling_rate = 0.5;
filter_half_width = 6; %default: 6;  % half width of filter (spatially and temporally)
subsampleRaw = 0;
%doMarginalize = 1; % default: 1     % option for specifying whether you want to capture just dynamics (1) or appearance and dynamics (0)
% G2H2OrG3 = 1;           % 0 = G2/H2, 1 = G3
% orientationSet = 2; % we used 2 a lot; 6, used in KTH and UCF sports experiment

%%

% pick the orientation set that we want to use (note: we can expan this to
% other orientation sets)
[nOrient, allOrientations, energyEps] = steerableFeatureSetProperties(orientationSet);


%% process each video
clear curVideoVol volG3 volSE tmpVol normVol
clear G2a_img G2b_img G2c_img G2d_img G2e_img G2f_img
clear H2a_img H2b_img H2c_img H2d_img H2e_img H2f_img H2g_img H2h_img H2i_img H2j_img

%     ref_filename = ref_videos{cVideo};
% 	dist_filename = dist_videos{cVideo};
% %  	v1 = Read_YUV([videoPath '/' ref_filename], frameSize);
%     fprintf('Reading in video %s\n', dist_filename);
%     tic
%  	videoVol =  Read_YUV([videoPath '/' dist_filename], frameSize);
%     frameNum = size(videoVol, 3);
%     fprintf('Done!'); toc
%     fprintf('Video length: %d frames\n', frameNum);
%     save videoVol videoVol
%     

%     load([videoPath '/' ref_filename(1:end-4) '.mat']);
%     ref_videoVol = videoVol(:,:, 1:30);
%     load([videoPath '/' dist_filename(1:end-4) '.mat']);
%     dist_videoVol = videoVol(:,:, 1:30);
    
    % file name for energies - should be changed/made automatic
%     storeFileName = [outputFolder, 'ETHPed2', ...
%                         '-eps-', num2str(epsilon),...
%                         '-dirSet-', num2str(orientationSet),...
%                         '-aggSize-', num2str(agg_size,'%02d'),...
%                         '-doMarg-', num2str(doMarginalize), ...
%                         '-G2G3-', num2str(G2H2OrG3), ...
%                         '.mat'];
           
   
    %% spatially subsample, if desired
    if subsampleRaw
        kernel = 1/16*[1 4 6 4 1]';
        curVideoVol = convn(convn(videoVol,...
            kernel, 'same'), kernel', 'same');
        curVideoVol = curVideoVol(1:2:end, 1:2:end, :);
    else
        curVideoVol = videoVol;
    end
    
    clear videoVol
    
   
    %% extract oriented energies with subsampling. 0 = G2/H2, 1 = G3
    if(G2H2OrG3 == 0)
        fprintf('Initializing the 3D G2 filters ...\n');
        tic;
        [G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img] = imgInit3DG2Sing(curVideoVol);
        fprintf('Done!'); toc
        
        fprintf('Initializing the 3D H2 filters ...\n');
        tic
        [H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img] = imgInit3DH2Sing(curVideoVol);
        fprintf('Done!'); toc
        
        tic
        fprintf('Computing the signatures with G2H2 ...\n');
        [volSE, volSE_n] = compute_signatures_with_G2H2_general(G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img, ...
                    H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img, ...
                    epsilon, allOrientations);              
        fprintf('Done!'); toc
        
        
    elseif(G2H2OrG3 == 1)
        fprintf('initializing 3D G3 ...\n');
        tic;
        volG3 = imgInit3DG3_with_options(curVideoVol,sampling_rate, filter_half_width);
        fprintf('Done!'); toc
        
        fprintf('computing signatures with G3 ...\n');
        tic
        [volSE, volSE_n] = compute_signatures_with_G3_general(volG3, epsilon, allOrientations, agg_size, doMarginalize);     
        fprintf('Done!'); toc;
        % volSE{1}   -- static channel
        % volSE{end} -- epsilon channel
        % sum the responses in the neighbourhood and normalized
    end   
%     fprintf('saving data to the disk ...\n')
%     tic
%     save(storeFileName, 'volSE', '-v7.3');     
%     fprintf('Done!'); toc
    


%%
% % can use the below for displaying the computed energies if everything is
% % still in memory
% currVidHists = single(cat(4, volSE{:}));
% theMin = min(currVidHists(:, :, :, 1:7));
% theMax = max(currVidHists(:, :, :, 1:7));
% theMin = min(theMin(:));
% theMax = max(theMax(:));
% 
%         currVidHists2 = single(cat(4, volSE{:}));
% theMin2 = min(currVidHists2(:, :, :, 1:7));
% theMax2 = max(currVidHists2(:, :, :, 1:7));
% theMin2 = min(theMin2(:));
% theMax2 = max(theMax2(:));
% 
% 
% for jj = 1:size(currVidHists2, 4)
%     for kk = 7:size(currVidHists2, 3)-7
%         %figure(1)
%         %imshow(currVidHists(:, :, kk, jj), [theMin theMax]);
%         
%         figure(2)
%         imshow(currVidHists2(:, :, kk, jj), [theMin2 theMax2]);
%         pause(0.1);
%     end
% end

%% play the original video
% videoSize = size(curVideoVol);
% temp = zeros(videoSize(1), videoSize(2), 1, videoSize(3));
% temp(:,:,1,:) = curVideoVol;
% temp = uint8(temp);
% mov = immovie(temp, colormap(gray));
% implay(mov)

% implay(curVideoVol);
