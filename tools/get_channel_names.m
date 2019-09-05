function selected_channel_names = get_channel_names(channel_names, index)


selected_channel_names = {channel_names{index}};
%return {channel_names{index}}

% remove 'ch_'
for i=1:length(selected_channel_names)
    loc_underscore = find(selected_channel_names{i}=='_');
    if ~isempty(loc_underscore)
        selected_channel_names{i} = selected_channel_names{i}(loc_underscore+1:end);
    end
end


