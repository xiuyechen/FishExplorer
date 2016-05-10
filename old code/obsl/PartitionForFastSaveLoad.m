
%% Partition test
nParts = 8;
t = zeros(1,nParts);
ix = round(linspace(1,size(CONST.CRZt,1),nParts+1));
for i = 1:nParts,
    tic;
    eval(['CRZt_' num2str(i) '= CONST.CRZt(ix(i):ix(i+1),:);']);
    save('CONST_F9_multi.mat',['CRZt_' num2str(i)'],'-v6','-append');
    t(i) = toc;
    i
    t(i)
end
sum(t)

%% reconstruct
tic

CRZt = zeros(size(CONST.CRZt));
num = 0;
for i = 1:8,
    eval(['len = size(CRZt_' num2str(i) ',1);']);
    eval(['CRZt(num+1:num+len,:) = CRZt_' num2str(i) ';']);
    num = num+len;
end

toc

%% for reals: partition and save
tic
const = CONST;
const = rmfield(const,'CRZt');
save('CONST_F9_fast.mat','const','-v6');

nParts = 8;
t = zeros(1,nParts);
ix = round(linspace(1,size(CONST.CRZt,1),nParts+1));
for i = 1:nParts,
    disp(num2str(i));
    eval(['CRZt_' num2str(i) '= CONST.CRZt(ix(i):ix(i+1),:);']);
    save('CONST_F9_fast.mat',['CRZt_' num2str(i)'],'-v6','-append');
end
toc

%% load
tic
load('CONST_F9_fast.mat','const');
load('CONST_F9_fast.mat','dimCR');
CRZt = zeros(dimCR);
num = 0;
for i = 1:8,
    load('CONST_F9_multi.mat',['CRZt_' num2str(i)']);
    eval(['len = size(CRZt_' num2str(i) ',1);']);
    eval(['CRZt(num+1:num+len,:) = CRZt_' num2str(i) ';']);
    num = num+len;
end
toc

%% reformat old mats
tic
data_dir = 'F:\Janelia2014';
for i_fish = 2:8,
    disp(['fish ' num2str(i_fish)]);
    fishdir = fullfile(data_dir,['CONST_F' num2str(i_fish) '.mat']);
    newfishdir = fullfile(data_dir,['CONST_F' num2str(i_fish) '_fast.mat']);
    load(fishdir,'CONST');
    const = CONST;
    const = rmfield(const,'CRZt');
    dimCR = size(CONST.CRZt);
    save(newfishdir,'const','-v6');
    save(newfishdir,'dimCR','-v6','-append');
    
    nParts = round(numel(CONST.CRZt)/42000000);
    t = zeros(1,nParts);
    ix = round(linspace(1,size(CONST.CRZt,1),nParts+1));
    for i = 1:nParts,
        disp(num2str(i));
        eval(['CRZt_' num2str(i) '= CONST.CRZt(ix(i):ix(i+1),:);']);
        save(newfishdir,['CRZt_' num2str(i)],'-v6','-append');
    end

end
toc


