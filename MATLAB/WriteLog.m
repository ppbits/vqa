%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writes a log 
%
% Name: Peng Peng
% Contact: dante.peng@gmail.com
% Date: Sept 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteLog(str)
timestamp = datestr(now);
fprintf('%s\t%s\n', timestamp, str);
fid = fopen('Logs.txt', 'at');
fprintf(fid, '%s\t%s\n', timestamp,str);
fclose(fid);
end
