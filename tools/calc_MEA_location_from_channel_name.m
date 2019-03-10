function xy = calc_MEA_location_from_channel_name(channel_name, ab, cd)

idx1 = str2num(channel_name(4));
idx2 = str2num(channel_name(5));


xy(2) = [idx1 1]*ab;
xy(1) = [idx2 1]*cd;


