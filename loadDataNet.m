
%% Choose set to load

% data in 20180828_GLM_CouplingFilter_SameRGCtype
base_folder_name = 'data/20180828_GLM_CouplingFilter_SameRGCtype';

% 10 Hz or 25 Hz data
StimInfoFileName = fullfile(base_folder_name, sprintf('StimInfo_8pix_163um_%dHz.mat',fps))
SpikeTrainFolderName = fullfile(base_folder_name, sprintf('SpikeTrain_20180828/%scell/%dHz/', CELL_TYPE, fps))


% % data in 20180820
% StimInfoFileName = '20180820/StimInfo_8pix_200um_10Hz.mat'
% SpikeTrainFolderName = sprintf('20180820/SpikeTrain_20180820/%scell/set%d/', CELL_TYPE, NUM_SET)


% % data in 20180724
% StimInfoFileName = '20180724/StimInfo_8pix_200um_10Hz.mat'
% SpikeTrainFolderName = sprintf('20180724/SpikeTrain_20180724/%scell/set%d/', CELL_TYPE, NUM_SET)

height=8,width=8

        
dtStim = 1 / fps;
dtSpike = 1 / fps;


%% load StimInfo and save to 'stim' of size (13 x 13 x 1800)
load(StimInfoFileName)

T = length(StimInfo);

stim = cell2mat(StimInfo);
%stim = reshape(stim, [T, height, width]);
% rearrange order
%stim = permute(stim, [2 3 1]);
size(stim)

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
    channelName = fileName(1:end-4);  % remove '.mat'
    if ~strcmp(channelName, 'A1a')
        channelNames{cnt} = channelName;    % save channel names (except A1a)
        cnt = cnt + 1;  % number of found channels
    end
    
    fprintf('Loading %s\n', fileName)
    load(fullfile(SpikeTrainFolderName, fileName))
    
    % save to txt file
    save([SpikeTrainFolderName channelName '.txt'], channelName, '-ascii')
    
end








                                 