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
    
    disp(sprintf('Loading %s', fileName))
    load([folderName fileName])
    
    % save to txt file
    save([folderName channelName '.txt'], channelName, '-ascii')
    
    %save(
%     if strcmp(fileName,'A1a.mat') % stim timing info
%         
% 
%         load([folderName fileName])
%     else
%         disp(sprintf('Loading %s', fileName))
%         load([folderName fileName])
%     end

end






                                 