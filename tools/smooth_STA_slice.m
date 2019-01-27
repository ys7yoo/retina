function STA_smoothed = smooth_STA_slice(STA, sig, height,width)

[T, num_pixels] = size(STA);

if nargin<2
    sig = 1;
end

if nargin<3
    height = sqrt(num_pixels);
    width = sqrt(num_pixels);
end



%% smooth each slice 
for t=1:T
    STA_smoothed(t,:,:) = imgaussfilt(reshape(STA(t,:),height,width), sig);
end

STA_smoothed = reshape(STA_smoothed, T, []);

