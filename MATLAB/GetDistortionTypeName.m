%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gets the distortion type name from the distortin type index
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Oct 17, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name = GetDistortionTypeName(database, index)
if(strcmp(database, 'LiveMobile'))
    if(index > 5 | index < -5)
        error('Distortion type index must be in the range of [-5, 5]');
    end

    names = { 
     'Compression',
     'Frame Freeze',
     'Rate Adaption',
     'Temporal Dynamics',
     'Wireless Channel Packet-loss'};

    if(index == 0)
        name = 'All';
    else if(index < 0)
        name = char(strcat('All (except for', {' '}, names{abs(index)}, ')'));
    else
        name = names{index};
    end
end
end
