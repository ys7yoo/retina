function [ab, cd] = fit_RF_center_onto_MEA(RFs, channel_names)

% extract found RFs
RF_center = [];
MEA_xy = [];
for n=1:length(RFs)
    if isfield(RFs{n}, 'mean')
        RF_center = [RF_center; RFs{n}.mean];
        
        x = str2num(channel_names{n}(4));
        y = str2num(channel_names{n}(5));
        MEA_xy = [MEA_xy; x, y];

    end
end



%% 
N = size(RF_center,1)

if N<2
    disp('not enough RFs')
    
    ab = nan;
    cd = nan;
    
    return
end

%% 
clf
plot(RF_center(:,1), RF_center(:,2), 'ok', 'markersize', 20)
hold on
for i=1:length(RF_center)
    text(RF_center(i,1), RF_center(i,2), sprintf('%d',MEA_xy(i,:)), 'HorizontalAlignment','center')
end


%%
clf
subplot(221)
plot(MEA_xy(:,1), RF_center(:,2),'o'); hold on
xlabel('channel index 1')
ylabel('RF_y')
subplot(222)
plot(MEA_xy(:,2), RF_center(:,1),'s')
xlabel('channel index 2')
ylabel('RF_x')


%% fit
ab = [MEA_xy(:,1) ones(N,1)]\RF_center(:,2)
cd = [MEA_xy(:,2) ones(N,1)]\RF_center(:,1)

subplot(221)
hold on 
plot([1 8], [[1 1]*ab, [8 1]*ab], 'r--')
title(sprintf(' y = %.1f x + %.1f', ab(1), ab(2)))
box off


subplot(222)
hold on 
plot([1 8], [[1 1]*cd, [8 1]*cd], 'r--')
title(sprintf(' y = %.1f x + %.1f', cd(1), cd(2)))
box off
%[8 1]*ab

% conversion matrix
A = [ab'; cd']

%% residual
subplot(223)
plot(MEA_xy(:,1), RF_center(:,2)-[MEA_xy(:,1) ones(N,1)]*ab, 'o')
xlabel('channel index 1')
ylabel('residual y')
box off

subplot(224)
plot(MEA_xy(:,2), RF_center(:,1)-[MEA_xy(:,2) ones(N,1)]*cd, 'o')
xlabel('channel index 2')
ylabel('residual x')
box off

return
