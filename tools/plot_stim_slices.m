function plot_stim_slices(stim, sta_num_samples, w, h)


if nargin<3
    dim = length(stim)/sta_num_samples;
    
    w = sqrt(dim);
    h = w;
end


stim_max = max(stim(:));
stim_min = min(stim(:));

st = reshape(stim, sta_num_samples, h, w);

c = 5;
r = ceil(sta_num_samples/c);





for i=1:sta_num_samples
    subplot(r,c,i)
    imagesc(reshape(st(i,:,:),h,w), [stim_min stim_max])
    %surf(reshape(st(i,:,:),h,w))
    colormap gray
    axis ij
end

return 


%% 
%sta_num_samples=10;
close all
idx = length(ev);
plot_stim_slices(U(:,idx), sta_num_samples)




% overlay RF
c = 5;
r = ceil(sta_num_samples/c);

n_RF = 1;
for i=1:sta_num_samples
    subplot(r,c,i); hold on
    switch RFs{n_RF}.type
        case 'ON'
            plot_ellipse(RFs{n_RF}.mean, RFs{n_RF}.cov, 'r-');
            tt=text(RFs{n_RF}.mean(1), RFs{n_RF}.mean(2), RFs{n_RF}.channel_name(4:end), 'HorizontalAlignment','center');
            tt.Color = [1 0 0];
        case 'OFF'
            plot_ellipse(RFs{n_RF}.mean, RFs{n_RF}.cov, 'b-');
            tt=text(RFs{n_RF}.mean(1), RFs{n_RF}.mean(2), RFs{n_RF}.channel_name(4:end), 'HorizontalAlignment','center')
            tt.Color = [0 0 1];
    end
end

%% 
[mask, XX, YY] = generate_ellipse_mask(RFs{n_RF}.mean, RFs{n_RF}.cov, 26, 26);
RF_area_in_pixel = sum(mask(:)>0)
mask_ext = repmat(mask(:)',sta_num_samples,1);
mask_ext = mask_ext(:);

plot_stim_slices(U(:,1).*mask_ext, sta_num_samples)


U_in_RF = reshape(U(mask_ext>0,1),10,[]);

xx = XX(mask>0);
yy = YY(mask>0);


%%

close all
plot(U_in_RF+repmat(0.05*(1:RF_area_in_pixel),sta_num_samples,1))



%% calc STC inside of RF


[sta_RF, stc_RF] = calc_STA_and_STC(stim(:,mask(:)>0), spike_train(:,n), sta_num_samples);
 
sta_RF_var = var(sta_RF,[],1);
[~, max_idx] = max(sta_RF_var)
[~, sorted_index] = sort(sta_RF_var, 'descend')

tic;
[u, d, ~] = svd(stc_RF);
toc



%Us = reshape(U(:,1),sta_num_samples,[]);
us = reshape(u(:,end),sta_num_samples,[]);
us_var = var(us, [], 1);

ev = diag(D);
    

close all
subplot(231)
plot(sta_RF)

subplot(232)
bar(sta_RF_var)

subplot(234)
plot(ev)
hold on
plot(1,ev(1), '*')
ylabel('eigen values')

subplot(235)
h = plot(us(:,sorted_index) - repmat(0.2*(1:RF_area_in_pixel),sta_num_samples,1) );
%h(max_idx).LineStyle = '--';
%plot(reshape(U(:,1),sta_num_samples,[]))
axis off

subplot(236)
%plot3(repmat((1:sta_num_samples)',1,25), repmat(xx',sta_num_samples,1), repmat(yy',sta_num_samples,1) + reshape(U(:,1),sta_num_samples,[])) 
bar(us_var)
xlabel('t')
ylabel('x')
zlabel('y')



%%
clf
loglog(sta_RF_var, us_var, '.'); xlabel('var(STA)');ylabel('var(STC)')


%% 
%plot_stim_slices(U(:,1), sta_num_samples)
plot_stim_slices(U(:,end), sta_num_samples)


