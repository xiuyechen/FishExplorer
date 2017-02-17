function file_masterdir = GetMasterFileDir
file_masterdir = 'C:\Users\Xiu\Dropbox (Personal)';
if ~exist(file_masterdir,'dir')
    file_masterdir = 'C:\Users\fish\Dropbox';
end

end