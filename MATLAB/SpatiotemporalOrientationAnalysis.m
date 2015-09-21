function [volSE, volSE_n] = SpatiotemporalOrientationAnalysis(...
    videoVol, allOrientations, G2H2OrG3, ...
    doMarginalize, filter_half_width)
if nargin < 5
    filter_half_width = 6;
end

if nargin < 4
    doMarginalize = 1;
end

if nargin < 3
    G2H2OrG3 = 1;           % G3
end


%% define necessary processing parameters
epsilon = 1000; %25 can adjust the eps value (from anywhere from very small (e.g., 1) to very larger (e.g., 10,000) depending on noise in video
padding = 50;       % not used now
agg_size = 5;% -1;      % for computing normalizing energy (5 is normally used)
sampling_rate = 0.5;
subsampleRaw = 0;
% filter_half_width = 6; %default: 6;  % half width of filter (spatially and temporally)
%doMarginalize = 1; % default: 1     % option for specifying whether you want to capture just dynamics (1) or appearance and dynamics (0)
% G2H2OrG3 = 1;           % 0 = G2/H2, 1 = G3
% orientationSet = 2; % we used 2 a lot; 6, used in KTH and UCF sports experiment




%% process each video
clear curVideoVol volG3 volSE tmpVol normVol
clear G2a_img G2b_img G2c_img G2d_img G2e_img G2f_img
clear H2a_img H2b_img H2c_img H2d_img H2e_img H2f_img H2g_img H2h_img H2i_img H2j_img

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
%     fprintf('initializing 3D G3 ...\n');
%     tic;
    volG3 = imgInit3DG3_with_options(curVideoVol,sampling_rate, filter_half_width);
%     fprintf('Done!'); toc
    
%     fprintf('computing signatures with G3 ...\n');
%     tic
    [volSE, volSE_n] = compute_signatures_with_G3_general(volG3, epsilon, allOrientations, agg_size, doMarginalize);
%     fprintf('Done!'); toc;
    % volSE{1}   -- static channel
    % volSE{end} -- epsilon channel
    % sum the responses in the neighbourhood and normalized
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G3   Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G3 = imgInit3DG3_with_options(vol, sampling_rate, filter_half_width)

if (nargin < 2)
    sampling_rate = 0.5;
    filter_half_width = 6;
elseif (nargin < 3)
    filter_half_width = 6;
end

t = sampling_rate*[-filter_half_width:filter_half_width];
C = 0.1840;

f1 = -4*C*t.*(-3 + 2*t.^2).*exp(-t.^2);
f2 = t.*exp(-t.^2);
f3 = -4*C*(-1 + 2*t.^2).*exp(-t.^2);
f4 = exp(-t.^2);
f5 = -8*C*t.*exp(-t.^2);

G3.a = aux_conv3d(vol, f1, f4, f4);
G3.b = aux_conv3d(vol, f3, f2, f4);
G3.c = aux_conv3d(vol, f2, f3, f4);
G3.d = aux_conv3d(vol, f4, f1, f4);
G3.e = aux_conv3d(vol, f3, f4, f2);
G3.f = aux_conv3d(vol, f5, f2, f2);
G3.g = aux_conv3d(vol, f4, f3, f2);
G3.h = aux_conv3d(vol, f2, f4, f3);
G3.i = aux_conv3d(vol, f4, f2, f3);
G3.j = aux_conv3d(vol, f4, f4, f1);

end


%% Called by imgInit3DG3_with_options
function out = aux_conv3d(vol, f1, f2, f3)
filter_length = length(f1);

out = convn(vol,        f1, 'same');   % x-direction
out = convn(out,        f2','same');   % y-direction
out = convn(out,reshape(f3,1,1,filter_length),'valid');   % z-direction

% out = imfilter(vol,        f1,'conv','same','symmetric');   % x-direction
% out = imfilter(out,        f2','conv','same','symmetric');   % y-direction
% out = imfilter(out,reshape(f3,1,1,filter_length),'conv','same','symmetric');   % z-direction
end



%% Function: compute G3 signatures
function [E, E_n] = compute_signatures_with_G3_general(G3, epsilon, orientations, agg_size, doMarginalize)

if ~exist('epsilon', 'var')
    epsilon = 500;
end

if ~exist('orientations', 'var')
    %if nargin < 3
    orientations = [0,  0, 1; % static channel is always the first
        1,  0, 1;
        -1,  0, 1;
        0,  1, 1;
        0, -1, 1;
        0,  1, 0];
end
nOrientations = size(orientations,1)+1;

if ~exist('agg_size', 'var') || isempty(agg_size)
    var_size = 13;
end

if doMarginalize <= 1
    nColumn = 1;
elseif doMarginalize == 2
    nColumn = 2;
end

E = cell(nOrientations, nColumn);

% compute base energies
for e = 1:nOrientations-1
    if(doMarginalize == 1)
        E{e, 1} = imgComputeSpacetimeOrientationEnergyG3(G3, orientations(e,1), orientations(e,2), orientations(e,3), agg_size);
    elseif(doMarginalize == 0)
        E{e, 1} = imgComputeSOEsG3_NoMarginalize(G3, orientations(e,:), agg_size);
    elseif doMarginalize == 2
        [E{e, 1}, E{e,2}]= imgComputeSpacetimeOrientationEnergyG3(G3, orientations(e,1), orientations(e,2), orientations(e,3), agg_size);
    end
end

E{nOrientations, 1} = epsilon*ones(size(E{1,1})); % last orientation is always epsilon
% last orientation is always epsilon
temp = cell(4,1);
for i = 1:4
    temp{i,1} = epsilon*ones(size(E{1,1}));
end
E{nOrientations, 2} = temp;

% normalize, if nt needed
E_n = cell(nOrientations, nColumn);
Esum = E{1, 1};
for e = 2:nOrientations
    Esum = Esum + E{e, 1};
end

for e = 1:nOrientations
    E_n{e, 1} = E{e, 1}./Esum;
end
for e = 1:nOrientations
    for sub_channel = 1:4 % four equally spaced directions on the plane
        E_n{e, 2}{sub_channel,1} = E{e,2}{sub_channel,1} ./ Esum;
    end
end
end


function [motion_energy] = imgComputeSOEsG3_NoMarginalize(G3, orientations, agg_size)

if (nargin < 3 || isempty(agg_size))
    agg_size = 13;
end

g = imgBinomialFilter(agg_size);

% if ~isempty(dir('imgSteer3DG3_fast.mex*'))
%     % use fast addition version
%     motion_energy = ...
%         imgSteer3DG3_fast(single(cos(0)*ua'     + sin(0)*ub'     ), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(1*pi/4)*ua'+ sin(1*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(2*pi/4)*ua'+ sin(2*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(3*pi/4)*ua'+ sin(3*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
% else
motion_energy = imgSteer3DG3(orientations, G3).^2;

%end

% perform aggregation only if kernel is specified
if agg_size > 0
    motion_energy = imfilter(imfilter(imfilter(motion_energy,g),g'),reshape(g,1,1,length(g)));
end

end


function [motion_energy, energy4channel] = imgComputeSpacetimeOrientationEnergyG3(G3, u, v, w, agg_size)
% N+1 equally spaced directions
energy4channel = cell(4,1);

if (nargin < 5 || isempty(agg_size))
    agg_size = 13;
end

g = imgBinomialFilter(agg_size);

n = [u; v; w];
n = n/norm(n);

e1 = [1; 0; 0];
e2 = [0; 1; 0];

if ( abs(acos(dot(n, e1)/norm(n))) > abs(acos(dot(n, e2)/norm(n))) )
    ua = cross(n,e1);
else
    ua = cross(n,e2);
end

ua = ua/norm(ua);

ub = cross(n,ua);
ub = ub/norm(ub);

% if you need more speed, should be able to uncomment this and it will use
% the mexx function (speeds things up roughly by a factor of 2)
if ~isempty(dir('../doesnotexist/imgSteer3DG3_fast.mex*')) % some time get Inf and NAN results
%     fprintf('fast');
    %    use fast addition version    
    energy4channel{1} = imgSteer3DG3_fast(single(cos(0)*ua'     + sin(0)*ub'     ), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    energy4channel{2} = imgSteer3DG3_fast(single(cos(1*pi/4)*ua'+ sin(1*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    energy4channel{3} = imgSteer3DG3_fast(single(cos(2*pi/4)*ua'+ sin(2*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    energy4channel{4} = imgSteer3DG3_fast(single(cos(3*pi/4)*ua'+ sin(3*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
    
%     motion_energy = ...
%         imgSteer3DG3_fast(single(cos(0)*ua'     + sin(0)*ub'     ), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(1*pi/4)*ua'+ sin(1*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(2*pi/4)*ua'+ sin(2*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2 + ...
%         imgSteer3DG3_fast(single(cos(3*pi/4)*ua'+ sin(3*pi/4)*ub'), G3.a,G3.b,G3.c,G3.d,G3.e,G3.f,G3.g,G3.h,G3.i,G3.j).^2;
else
%     fprintf('no fast');
    energy4channel{1} =  imgSteer3DG3(cos(0)*ua'      + sin(0)*ub'     , G3).^2 ;
    energy4channel{2} =   imgSteer3DG3(cos(1*pi/4)*ua' + sin(1*pi/4)*ub', G3).^2;
    energy4channel{3} =  imgSteer3DG3(cos(2*pi/4)*ua' + sin(2*pi/4)*ub', G3).^2;
    energy4channel{4} =   imgSteer3DG3(cos(3*pi/4)*ua' + sin(3*pi/4)*ub', G3).^2;
    
%         motion_energy = imgSteer3DG3(cos(0)*ua'      + sin(0)*ub'     , G3).^2 + ...
%         imgSteer3DG3(cos(1*pi/4)*ua' + sin(1*pi/4)*ub', G3).^2 + ...
%         imgSteer3DG3(cos(2*pi/4)*ua' + sin(2*pi/4)*ub', G3).^2 + ...
%         imgSteer3DG3(cos(3*pi/4)*ua' + sin(3*pi/4)*ub', G3).^2;
end

motion_energy = energy4channel{1} + energy4channel{2} + energy4channel{3} + energy4channel{4};

% perform aggregation only if kernel is specified
if agg_size > 0
    motion_energy = imfilter(imfilter(imfilter(motion_energy,g),g'),reshape(g,1,1,length(g)));
    for i = 1:4 
        energy4channel{1}= imfilter(imfilter(imfilter(energy4channel{1},g),g'),reshape(g,1,1,length(g)));
    end
end

end


function [vol] = imgSteer3DG3(direction, G3)

a = direction(1);
b = direction(2);
c = direction(3);

vol =     (a^3)*G3.a ...
    + 3*(a^2)*b*G3.b ...
    + 3*a*(b^2)*G3.c ...
    +     (b^3)*G3.d ...
    + 3*(a^2)*c*G3.e ...
    +   6*a*b*c*G3.f ...
    + 3*(b^2)*c*G3.g ...
    + 3*a*(c^2)*G3.h ...
    + 3*b*(c^2)*G3.i ...
    +     (c^3)*G3.j;
end


% Description: Generates a normalized binomial filter where n represents
% the desired number of taps
function bfilter = imgBinomialFilter(n)

n = n - 1;

bfilter = [];
for i = 0:n
    bfilter = [bfilter nchoosek(n,i)];
end

% --- normalization ---
bfilter = bfilter/sum(bfilter);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% G2H2 Filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E, E_n] = compute_signatures_with_G2H2_general(G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img, ...
    H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img, ...
    epsilon, orientations)

if ~exist('epsilon', 'var')
    epsilon = 500;
end

% add one more orientation for epsilon channel
nOrientations = size(orientations,1)+1;

E = cell(nOrientations, 1);

% compute base energies
for e = 1:nOrientations-1
    E{e} = imgSteer3DG2(orientations(e, :), G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img).^2 + ...
        imgSteer3DH2(orientations(e, :), H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img).^2;
end
E{nOrientations} = epsilon*ones(size(E{1})); % last orientation is always epsilon

% sum across channel to form the denominator
Esum = E{1};
for e = 2:nOrientations
    Esum = Esum + E{e};
end

% normalize by the sum of the consort
E_n = cell(nOrientations, 1);
for e = 1:nOrientations
    E_n{e} = E{e}./Esum;
end

end


function [G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img] = imgInit3DG2(img)
% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---

SAMPLING_RATE = 0.67;
N = (2/sqrt(3))*(2/pi)^(3/4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms t;
% f1 = N*(2*t^2 - 1)*exp(-t^2);
% f2 = exp(-t^2);
% f3 = 2*N*t*exp(-t^2);
% f4 = t*exp(-t^2);
%
% i = SAMPLING_RATE*[-4:4];
% filter_size = length(i);
%
% basis = [subs(f1,t,i); subs(f2,t,i); subs(f3,t,i); subs(f4,t,i)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KC modification to get rid of syms and subs (much faster and works on
% Laplace)
i = SAMPLING_RATE*[-4:4];
filter_size = length(i);

t = i;

f1 = N*(2*t.^2 - 1).*exp(-t.^2);
f2 = exp(-t.^2);
f3 = 2*N*t.*exp(-t.^2);
f4 = t.*exp(-t.^2);
basis = [f1; f2; f3; f4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G2a_img = single(convn(img, reshape(basis(1,:),1,filter_size),'same'));   % x-direction
G2a_img = single(convn(G2a_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
G2a_img = single(convn(G2a_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

G2b_img = single(convn(img, reshape(basis(3,:),1,filter_size),'same'));   % x-direction
G2b_img = single(convn(G2b_img,reshape(basis(4,:),filter_size,1),'same'));   % y-direction
G2b_img = single(convn(G2b_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

G2c_img = single(convn(img, reshape(basis(2,:),1,filter_size),'same'));   % x-direction
G2c_img = single(convn(G2c_img,reshape(basis(1,:),filter_size,1),'same'));   % y-direction
G2c_img = single(convn(G2c_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

G2d_img = single(convn(img, reshape(basis(3,:),1,filter_size),'same'));   % x-direction
G2d_img = single(convn(G2d_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
G2d_img = single(convn(G2d_img,reshape(basis(4,:),1,1,filter_size),'same')); % z-direction

G2e_img = single(convn(img, reshape(basis(2,:),1,filter_size),'same'));   % x-direction
G2e_img = single(convn(G2e_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
G2e_img = single(convn(G2e_img,reshape(basis(4,:),1,1,filter_size),'same')); % z-direction

G2f_img = single(convn(img, reshape(basis(2,:),1,filter_size),'same'));   % x-direction
G2f_img = single(convn(G2f_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
G2f_img = single(convn(G2f_img,reshape(basis(1,:),1,1,filter_size),'same')); % z-direction
end




function [H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img] = imgInit3DH2Sing(img)
% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---

SAMPLING_RATE = 0.67;
N = 0.877776;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms t;
%
% f1 = N*(t^3 - 2.254*t)*exp(-t^2);
% f2 = N*(t^2 - 0.751333)*exp(-t^2);
% f3 = exp(-t^2);
% f4 = N*t*exp(-t^2);
% f5 = t*exp(-t^2);
%
% i = SAMPLING_RATE*[-4:4];
% filter_size = length(i);
%
% basis = [subs(f1,t,i); subs(f2,t,i); subs(f3,t,i); subs(f4,t,i); subs(f5,t,i)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KC modification to get rid of syms and subs
i = SAMPLING_RATE*[-4:4];
filter_size = length(i);

t = i;

f1 = N*(t.^3 - 2.254.*t).*exp(-t.^2);
f2 = N*(t.^2 - 0.751333).*exp(-t.^2);
f3 = exp(-t.^2);
f4 = N*t.*exp(-t.^2);
f5 = t.*exp(-t.^2);

basis = [f1; f2; f3; f4; f5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H2a_img = single(convn(img,    reshape(basis(1,:),1,filter_size),'same'));   % x-direction
H2a_img = single(convn(H2a_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2a_img = single(convn(H2a_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2b_img = single(convn(img,    reshape(basis(2,:),1,filter_size),'same'));   % x-direction
H2b_img = single(convn(H2b_img,reshape(basis(5,:),filter_size,1),'same'));   % y-direction
H2b_img = single(convn(H2b_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2c_img = single(convn(img,    reshape(basis(5,:),1,filter_size),'same'));   % x-direction
H2c_img = single(convn(H2c_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
H2c_img = single(convn(H2c_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2d_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2d_img = single(convn(H2d_img,reshape(basis(1,:),filter_size,1),'same'));   % y-direction
H2d_img = single(convn(H2d_img,reshape(basis(3,:),1,1,filter_size),'same')); % z-direction

H2e_img = single(convn(img,    reshape(basis(2,:),1,filter_size),'same'));   % x-direction
H2e_img = single(convn(H2e_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2e_img = single(convn(H2e_img,reshape(basis(5,:),1,1,filter_size),'same')); % z-direction

H2f_img = single(convn(img,    reshape(basis(4,:),1,filter_size),'same'));   % x-direction
H2f_img = single(convn(H2f_img,reshape(basis(5,:),filter_size,1),'same'));   % y-direction
H2f_img = single(convn(H2f_img,reshape(basis(5,:),1,1,filter_size),'same')); % z-direction

H2g_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2g_img = single(convn(H2g_img,reshape(basis(2,:),filter_size,1),'same'));   % y-direction
H2g_img = single(convn(H2g_img,reshape(basis(5,:),1,1,filter_size),'same')); % z-direction

H2h_img = single(convn(img,    reshape(basis(5,:),1,filter_size),'same'));   % x-direction
H2h_img = single(convn(H2h_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2h_img = single(convn(H2h_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

H2i_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2i_img = single(convn(H2i_img,reshape(basis(5,:),filter_size,1),'same'));   % y-direction
H2i_img = single(convn(H2i_img,reshape(basis(2,:),1,1,filter_size),'same')); % z-direction

H2j_img = single(convn(img,    reshape(basis(3,:),1,filter_size),'same'));   % x-direction
H2j_img = single(convn(H2j_img,reshape(basis(3,:),filter_size,1),'same'));   % y-direction
H2j_img = single(convn(H2j_img,reshape(basis(1,:),1,1,filter_size),'same')); % z-direction
end



function [img] = imgSteer3DG2(direction, G2a_img, G2b_img, G2c_img, G2d_img, G2e_img, G2f_img)

% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---

a = direction(1);
b = direction(2);
c = direction(3);

img = (a^2)*G2a_img ...
    + 2*a*b*G2b_img ...
    + (b^2)*G2c_img ...
    + 2*a*c*G2d_img ...
    + 2*b*c*G2e_img ...
    + (c^2)*G2f_img;

end



function [img] = imgSteer3DH2(direction, H2a_img, H2b_img, H2c_img, H2d_img, H2e_img, H2f_img, H2g_img, H2h_img, H2i_img, H2j_img)

% ---
% Written by: Konstantinos G. Derpanis
% PhD Candidate
% York University
% Toronto, Ontario, Canada
% Bugs/Questions email me at: kosta@cs.yorku.ca
% see York Universtiy Technical Report CS-2004-05 for analytic details
%   located here: http://www.cs.yorku.ca/techreports/2004/CS-2004-05.html
% ---
a = direction(1);
b = direction(2);
c = direction(3);

img = (a^3)*H2a_img ...
    + 3*(a^2)*b*H2b_img ...
    + 3*a*(b^2)*H2c_img ...
    + (b^3)*H2d_img ...
    + 3*(a^2)*c*H2e_img ...
    + 6*a*b*c*H2f_img ...
    + 3*(b^2)*c*H2g_img ...
    + 3*a*(c^2)*H2h_img ...
    + 3*b*(c^2)*H2i_img ...
    + (c^3)*H2j_img;
end






