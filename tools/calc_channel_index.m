function index = calc_channel_index(channel_names, selected_channel_names)

if ~iscell(selected_channel_names) % single string
    selected_channel_names={selected_channel_names};
end
    
index = [];
for n = 1:length(selected_channel_names)
    channel_name_to_query = ['ch_' selected_channel_names{n}];
    idx = find(ismember(channel_names, channel_name_to_query));
    if isempty(idx)
        disp(sprintf('cannot find channel %s',channel_name_to_query))
    else
        
        index = [index idx];
    end
end

