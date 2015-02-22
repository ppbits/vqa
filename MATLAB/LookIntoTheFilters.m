% LookIntoTheFilters
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

allFilters = zeros(size(allOrient,1)*4, size(allOrient,2));
for theta = 1:size(allOrient,1)
    n = allOrient(theta,:)';
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
    
    n_k1 = cos(0)*ua'      + sin(0)*ub';
    n_k2 = cos(1*pi/4)*ua' + sin(1*pi/4)*ub';
    n_k3 = cos(2*pi/4)*ua' + sin(2*pi/4)*ub';
    n_k4 = cos(3*pi/4)*ua' + sin(3*pi/4)*ub';
    
    allFilters(1+(theta-1)*4: theta*4, :) = [n_k1; n_k2; n_k3,; n_k4];
end

[~, ind] = sort(allFilters(:,1));
allFilters_sorted = allFilters(ind,:);