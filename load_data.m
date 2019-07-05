function [stim, spike_train, channel_names, exp_param] = load_data(base_folder_name)

%% parse info from StimInfo file
find_stim_info_filename = dir(fullfile(base_folder_name,'Stiminfo*'));

num_stim_info_files=length(find_stim_info_filename);
if num_stim_info_files>1
    % there may be multiple files!
    disp('There are multiple files in the folder.')
    for i=1:num_stim_info_files
        disp(sprintf("%d: %s", i, find_stim_info_filename(i).name))
    end
    idx = input('Choose the index from above. ');
    stim_info_filename = find_stim_info_filename(idx).name;
else
    stim_info_filename = find_stim_info_filename.name;
    
end


%% parse info from StimTrain folder
spike_train_folder_name = dir(fullfile(base_folder_name,'SpikeTrain_*'));
spike_train_folder_name=spike_train_folder_name.name;


[num_pixels_per_dim, pixel_size, sampling_rate] = parse_stim_info_filename(stim_info_filename);

exp_param.num_pixels_per_dim = num_pixels_per_dim;
exp_param.pixel_size = pixel_size;
exp_param.sampling_rate = sampling_rate;

%% load 'StimInfo'
load(fullfile(base_folder_name,stim_info_filename));

% T = length(StimInfo);

stim = cell2mat(StimInfo);
%size(stim) % 9000 x num_pixels_per_dim x num_pixels_per_dim


%%  plot some sample stimulus
%clf
figure
subplot(141)
imshow(reshape(stim(1,:),[exp_param.num_pixels_per_dim, exp_param.num_pixels_per_dim]))
%xlabel('x');ylabel('y');
axis xy

subplot(142)
imshow(reshape(stim(2,:),[exp_param.num_pixels_per_dim, exp_param.num_pixels_per_dim]))
%xlabel('x');ylabel('y');
axis xy

subplot(143)
imshow(reshape(stim(3,:),[exp_param.num_pixels_per_dim, exp_param.num_pixels_per_dim]))
%xlabel('x');ylabel('y');
axis xy

subplot(144)
imshow(reshape(stim(4,:),[exp_param.num_pixels_per_dim, exp_param.num_pixels_per_dim]))
%xlabel('x');ylabel('y');
axis xy

%%  load SpikeTrain Dataset & save txt
exp_param.spike_train_folder_name = fullfile(base_folder_name, sprintf('%s/SpikeTrain_ND2_%dpix_%dum_%dHz', spike_train_folder_name, exp_param.num_pixels_per_dim, exp_param.pixel_size, exp_param.sampling_rate));
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
    for i = 2:length(A1a)
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

figure
%clf
imagesc(spike_train', [0 1]); colormap gray
xlabel('t')
ylabel('channel')