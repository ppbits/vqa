%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nOrient, allOrient, energyEps] = SteerableFeatureSetProperties(featureSet)
% return various properties accosiated with the featureSet used for
% steerable filter responses histogram

energyEps = 1e5;

r = (sqrt(5)+1)/2;            
% choose the 10 feature set
if(featureSet == 0)
    allOrient = [0 r 1/r; ...                         
        -1/r 0 r; ...               
        r 1/r 0; ...               
        0, -r, 1/r; ...             
        -1/r, 0, -r; ...            
        -r, 1/r, 0; ...            
        1, 1, -1; ...               
        -1, 1, 1; ...               
        -1, 1, -1; ...              
        -1, -1, -1];                
    % normalize vectors
    allOrient = allOrient./repmat(sqrt(sum(allOrient.^2, 2)), [1 3]);

% choose the 36 feature sset
elseif(featureSet == 1)

    allOrient = [ ...
        % XT first
        1, 0, 0; ...                    % vertical structures in the XT domain (static)
        -1/sqrt(2), 0, 1/sqrt(2); ...   % rightward motion (diagonal structure)
        0, 0, 1; ...                    % horizontal structures in the XT domain (flicker)
        1/sqrt(2), 0, 1/sqrt(2); ...    % left motion (diagonal structure)

        % now do YT
        0, 1, 0; ...                    % vertical structures in the YT domain (static)
        0, 1/sqrt(2), -1/sqrt(2); ...   % down motion (diagonal structure)
         0, 0, 1; ...                    % horizontal structures in the YT domain (flicker)
        0, 1/sqrt(2), 1/sqrt(2); ...    % up motion (diagonal structure)

        % now do XY
        1, 0, 0; ...                    % vertical structures in the XY domain       
        -1/sqrt(2), 1/sqrt(2), 0; ...   % diagonal structure in XY
        0, 1, 0; ...                    % horizontal structures in the XY domain                    
        1/sqrt(2), 1/sqrt(2), 0; ...    % diagonal structure in XY                    
        ];   

    
    %% The features below indicate the normals of the planes in the frequency domain 
    % 
% 13-direction non-repetitive feature set (Mikhail)
elseif(featureSet == 2)
    allOrient = [...
        1, 0, 0; ...    %flicker horizontal     1
        1, 1, 0; ...    %flicker diagonal a     2
        0, 1, 0; ...    %flicker vertical       3 
        -1, 1, 0; ...   %flicker diagonal b     4
        1, 1, 1; ...    %down right             5    
        0, 1, 1; ...    %up                     6       
        -1, 1, 1; ...   %down left              7 
        1, 0, 1; ...    %right                  8            
        0, 0, 1; ...    %static                 9   
        -1, 0, 1; ...   %left                   10   
        1, -1, 1; ...   %up right               11
        0, -1, 1; ...   %down                   12   
        -1, -1, 1; ...  %up left                13     
        ];
    % normalize vectors
    allOrient = allOrient./repmat(sqrt(sum(allOrient.^2, 2)), [1 3]);
    
elseif(featureSet == 3)
    allOrient =[0 0 1                           %static                 1
                1 0 0                           %flicker horizontal     2
                0 1 0                           %flicker vertical       3
                1 1 0                           %flicker diagonal a     4
                -1 1 0                          %flicker diagonal b     5
                0.5 0 1                         %right slow             6
                1 0 1                           %right                  7
                2 0 1                           %right fast             8
                -0.5 0 1                        %left slow              9
                -1 0 1                          %left                   10
                -2 0 1                          %left fast              11
                0 0.5 1                         %up slow                12
                0 1 1                           %up                     13
                0 2 1                           %up fast                14
                0 -0.5 1                        %down slow              15
                0 -1 1                          %down                   16
                0 -2 1                          %down fast              17
                0.5/sqrt(2) 0.5/sqrt(2) 1       %right-down slow        18
                1/sqrt(2) 1/sqrt(2) 1           %right-down             19
                2/sqrt(2) 2/sqrt(2) 1           %right-down fast        20
                0.5/sqrt(2) -0.5/sqrt(2) 1      %right-up slow          21
                1/sqrt(2) -1/sqrt(2) 1          %right-up               22
                2/sqrt(2) -2/sqrt(2) 1          %right-up fast          23
                -0.5/sqrt(2) 0.5/sqrt(2) 1      %left-down slow         24
                -1/sqrt(2) 1/sqrt(2) 1          %left-down              25
                -2/sqrt(2) 2/sqrt(2) 1          %left-down fast         26
                -0.5/sqrt(2) -0.5/sqrt(2) 1     %left-up slow           27 
                -1/sqrt(2) -1/sqrt(2) 1         %left-up                28
                -2/sqrt(2) -2/sqrt(2) 1];       %left-up fast           29
            
elseif(featureSet == 4)
    allOrient = [0,  0, 1; % static channel is always first     1
                 1,  0, 1;                      %right          2
                -1,  0, 1;                      %left           3
                 0,  1, 1;                      % down          4
                 0, -1, 1;                      % up            5             
                 0,  1, 0];                     % flicker       6
elseif(featureSet == 5)
    allOrient =[0 0 1                           %static                 1
                1 0 0                           %flicker horizontal     2
                0 1 0                           %flicker vertical       3                
                1 0 1                           %right                  4
                3 0 1                           %right fast             5               
                -1 0 1                          %left                   6
                -3 0 1                          %left fast              7
                0 1 1                           %up                     8
                0 3 1                           %up fast                9
                0 -1 1                          %down                   10
                0 -3 1                          %down fast              11
                ];
elseif(featureSet == 6)
    allOrient =[0 0 1                           %static                 1
                1 0 0                           %flicker horizontal     2
                0 1 0                           %flicker vertical       3                
                1 0 1                           %right                  4
                3 0 1                           %right fast             5               
                -1 0 1                          %left                   6
                -3 0 1                          %left fast              7
                0 1 1                           %up                     8
                0 3 1                           %up fast                9
                0 -1 1                          %down                   10
                0 -3 1                          %down fast              11
                1 1 0                           %flicker diagonal       12
                1 -1 0                          %flicker diagonal other 13                
                1 1 1                           %down right             14
                3 3 1                           %down right fast        15               
                -1 1 1                          %down left              16
                -3 3 1                          %down left fast         17
                1 -1 1                          %up right               18
                3 -3 1                          %up right fast          19
                -1 -1 1                         %up left                21
                -3 -3 1                         %up left fast           22
                ];
% similar to 4, but tuned to fast motions
elseif(featureSet == 7)
    allOrient = [0,  0, 1; % static channel is always first     1
                 3,  0, 1;                      % fast right    2
                -3,  0, 1;                      % fast left     3
                 0,  3, 1;                      % fast down     4
                 0, -3, 1;                      % fast up       5             
                 0,  1, 0];                     % flicker       6
end   




nOrient = size(allOrient, 1);

% make all vectors point in the positive t-direction
for cOrient = 1:nOrient
    if allOrient(cOrient,3) < 0
        allOrient(cOrient,:) = -allOrient(cOrient,:);
    end
end
    
