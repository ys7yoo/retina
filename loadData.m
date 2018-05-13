%% load Dataset & save txt
clear
folderName = 'SpikeTrain_20180404/SpikeTrain_Num1_13pix_100um_10Hz/'
%folderName = 'SpikeTrain_20180404/SpikeTrain_Num2_26pix_50um_10Hz/'
%folderName = 'SpikeTrain_20180404/SpikeTrain_Num3_13pix_100um_20Hz/'
%folderName = 'SpikeTrain_20180404/SpikeTrain_Num4_26pix_50um_20Hz/'

matFiles = dir([folderName '*.mat']);
matFiles = {matFiles.name};

N = length(matFiles);

for i = 1:N
    
    fileName = matFiles{i};    
    channelName = fileName(1:end-4);
    
    fprintf('Loading %s\n', fileName)
    load([folderName fileName])
    
    % save to txt file
    save([folderName channelName '.txt'], channelName, '-ascii')
    
end






                                 