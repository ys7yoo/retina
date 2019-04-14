function [stim, spike_train, channel_names, exp_param] = load_data(folder_name_by_date)

%% 
exp_param.folder_name_by_date = folder_name_by_date;
base_folder_name = fullfile('coupling_data', folder_name_by_date);

%% parse info from StimInfo file
find_stim_info_filename = dir(fullfile(base_folder_name,'Stiminfo*'));
stim_info_filename = find_stim_info_filename.name;


[num_pixels, distance_between_electrodes, sampling_rate] = parse_stim_info_filename(stim_info_filename);

exp_param.num_pixels = num_pixels;
exp_param.distance_between_electrodes = distance_between_electrodes;
exp_param.sampling_rate = sampling_rate;

%% load 'StimInfo'
load(fullfile(base_folder_name,stim_info_filename));

% T = length(StimInfo);

stim = cell2mat(StimInfo);

%size(stim) % 9000 x num_pixels x num_pixels


%%  plot some sample stimulus
clf
subplot(131)
imshow(reshape(stim(1,:),[exp_param.num_pixels, exp_param.num_pixels]))
xlabel('x');ylabel('y');axis xy

subplot(132)
imshow(reshape(stim(2,:),[exp_param.num_pixels, exp_param.num_pixels]))
xlabel('x');ylabel('y');axis xy

subplot(133)
imshow(reshape(stim(3,:),[exp_param.num_pixels, exp_param.num_pixels]))
xlabel('x');ylabel('y');axis xy

%%  load SpikeTrain Dataset & save txt
exp_param.spike_train_folder_name = fullfile(base_folder_name, sprintf('SpikeTrain_%s/SpikeTrain_ND2_%dpix_%dum_%dHz', exp_param.folder_name_by_date, exp_param.num_pixels, exp_param.distance_between_electrodes, exp_param.sampling_rate));
exp_param.spike_train_folder_name = [exp_param.spike_train_folder_name '/']

matFiles = dir([exp_param.spike_train_folder_name '*.mat']);
if isempty(matFiles)
    error(sprintf('Cannot load spike train data from %s', exp_param.spike_train_folder_name))
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
    load(fullfile(exp_param.spike_train_folder_name, fileName))
    
    % save to txt file
    save([exp_param.spike_train_folder_name channel_name '.txt'], channel_name, '-ascii')
    
end


%% For further analysis, need to convert spike time to spike train

N = length(channel_names); % number of neurons 

binStim = size(stim,1);
assert (binStim==length(A1a))
spike_train = zeros(binStim,N);
%spike_train = sparse(binStim,N);

for n = 1:N
    
    channel_name = channel_names{n};
    disp(channel_name)
    spike_time=eval(channel_name);
    
    
    % let's check for each stim time bin
    for i = 2:binStim
        t0 = A1a(i-1);
        t1 = A1a(i);

        % find stim time that occured during t0 and t1
        idx = find(spike_time>t0 & spike_time<=t1);

        spike_train(i,n) = length(idx);  % multiple spikes may occur in a bin
%         if ~isempty(idx)
%             %disp('found')
%             spike_train(i,n) = 1;
%         end

    end
end


clf
imagesc(spike_train', [0 1]); colormap gray


