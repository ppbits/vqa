% Write Logs to screen and file.
%
% Contact: dante.peng@gmail.com
% Copyright (c) Peng Peng
function WriteLog(str)
timestamp = datestr(now);
fprintf('%s\t%s\n', timestamp, str);
fid = fopen('Logs.txt', 'at');
fprintf(fid, '%s\t%s\n', timestamp,str);
fclose(fid);
end
