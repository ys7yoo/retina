function [num_pixels, pixel_size, sampling_rate] = parse_stim_info_filename(stim_info_filename)

strs = strsplit(stim_info_filename(1:end-4),'_');

num_pixels = str2num(strs{2}(1:end-3));
pixel_size = str2num(strs{3}(1:end-2));
sampling_rate = str2num(strs{4}(1:end-2));