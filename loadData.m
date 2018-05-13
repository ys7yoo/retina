
%% Choose experiment to load
clear

NUM_EXP = 1

switch NUM_EXP
    case 1        
        StimInfoFileName = '20180404/StimInfo_Num1_13pix_100um_10Hz.mat'
        SpikeTrainFolderName = '20180404/SpikeTrain_20180404/SpikeTrain_Num1_13pix_100um_10Hz/'
        
        height=13,width=13
        fps = 10 %  set manually according to the SpikeTrainFolderName        
    case 2
        StimInfoFileName = '20180404/StimInfo_Num2_26pix_50um_10Hz.mat'
        SpikeTrainFolderName = '20180404/SpikeTrain_20180404/SpikeTrain_Num2_26pix_50um_10Hz/'        
        
        height=26,width=26
        fps = 10 %  set manually according to the SpikeTrainFolderName        
    case 3
        StimInfoFileName = '20180404/StimInfo_Num3_13pix_100um_20Hz.mat'
        SpikeTrainFolderName = '20180404/SpikeTrain_20180404/SpikeTrain_Num3_13pix_100um_20Hz/'
       
        height=13,width=13
        fps = 20 %  set manually according to the SpikeTrainFolderName        
    case 4
        StimInfoFileName = '20180404/StimInfo_Num4_26pix_50um_20Hz.mat'
        SpikeTrainFolderName = '20180404/SpikeTrain_20180404/SpikeTrain_Num4_26pix_50um_20Hz/'
                
        height=26,width=26
        fps = 20 %  set manually according to the SpikeTrainFolderName        
        
    otherwise
        error('Cannot load data')
end


%% load StimInfo and save to 'stim' of size (13 x 13 x 1800)
load(StimInfoFileName)

T = length(StimInfo);

stim = cell2mat(StimInfo);
stim = reshape(stim, [T, height, width]);
% rearrange order
stim = permute(stim, [2 3 1]);
size(stim)

% plot some sample stimulus
clf
subplot(131)
imshow(stim(:,:,1))
subplot(132)
imshow(stim(:,:,2))
subplot(133)
imshow(stim(:,:,3))



%%  load SpikeTrain Dataset & save txt

matFiles = dir([SpikeTrainFolderName '*.mat']);
matFiles = {matFiles.name};

for i = 1:length(matFiles)
    
    fileName = matFiles{i};    
    channelName = fileName(1:end-4);
    channelNames{i} = channelName;    % save channel names 
    
    fprintf('Loading %s\n', fileName)
    load([SpikeTrainFolderName fileName])
    
    % save to txt file
    save([SpikeTrainFolderName channelName '.txt'], channelName, '-ascii')
    
end








                                 