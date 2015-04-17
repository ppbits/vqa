% Play ground
function WriteLog(str)
fprintf('%s\n', str);
fid = fopen('Logs.txt', 'at');
fprintf(fid, '%s\t%s\n', datestr(now), str);
fclose(fid);
end