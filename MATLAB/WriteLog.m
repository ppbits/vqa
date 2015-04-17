% Write Logs to screen and file.
%
% Contact: dante.peng@gmail.com
% Copyright (c) Peng Peng
function WriteLog(str)
fprintf('%s\n', str);
fid = fopen('Logs.txt', 'at');
fprintf(fid, '%s\t%s\n', datestr(now), str);
fclose(fid);
end