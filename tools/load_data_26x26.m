

%% load 26x26 stim results from 20180905_26x26

% data in 20180828_GLM_CouplingFilter_SameRGCtype
base_folder_name = 'data/20180905_26x26';

fps = 30
StimInfoFileName = fullfile(base_folder_name, sprintf('StimInfo_26pix_50um_%dHz.mat',fps))
SpikeTrainFolderName = fullfile(base_folder_name, sprintf('SpikeTrain_20180905/SpikeTrain_ND2_26pix_50um_%dHz/', fps))


% % data in 20180820
% StimInfoFileName = '20180820/StimInfo_8pix_200um_10Hz.mat'
% SpikeTrainFolderName = sprintf('20180820/SpikeTrain_20180820/%scell/set%d/', CELL_TYPE, NUM_SET)


% % data in 20180724
% StimInfoFileName = '20180724/StimInfo_8pix_200um_10Hz.mat'
% SpikeTrainFolderName = sprintf('20180724/SpikeTrain_20180724/%scell/set%d/', CELL_TYPE, NUM_SET)

height=26,width=26

        
dtStim = 1 / fps;
dtSpike = 1 / fps;


%% load StimInfo and save to 'stim' of size (13 x 13 x 1800)
load(StimInfoFileName)

T = length(StimInfo);

stim = cell2mat(StimInfo);
%stim = reshape(stim, [T, height, width]);
% rearrange order
%stim = permute(stim, [2 3 1]);
size(stim)                       % 27000         676


% plot some sample stimulus
clf
subplot(131)
imshow(reshape(stim(1,:),[height, width]))
xlabel('x');ylabel('y');axis xy

subplot(132)
imshow(reshape(stim(2,:),[height, width]))
xlabel('x');ylabel('y');axis xy

subplot(133)
imshow(reshape(stim(3,:),[height, width]))
xlabel('x');ylabel('y');axis xy


%%  load SpikeTrain Dataset & save txt

matFiles = dir([SpikeTrainFolderName '*.mat']);
if isempty(matFiles)
    error(sprintf('Cannot load spike train data from %s',SpikeTrainFolderName))
end
matFiles = {matFiles.name};

cnt = 1;
for i = 1:length(matFiles)
    
    fileName = matFiles{i};    
    channel_name = fileName(1:end-4);  % remove '.mat'
    if ~strcmp(channel_name, 'A1a')
        channel_names{cnt} = channel_name;    % save channel names (except A1a)
        cnt = cnt + 1;  % number of found channels
    end
    
    fprintf('Loading %s\n', fileName)
    load(fullfile(SpikeTrainFolderName, fileName))
    
    % save to txt file
    save([SpikeTrainFolderName channel_name '.txt'], channel_name, '-ascii')
    
end








                                 