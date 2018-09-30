%% Fit GLM with coupling filters!

%% set path to GLM tools

% This is where the GLMspiketools is installed
% You can get the package from https://github.com/ys7yoo/GLMspiketools
basedir = [getenv('HOME') '/src/GLMspiketools']

% Add a bunch sub-directories (with absoluate path names)
addpath([basedir '/glmtools_fitting/']);
addpath([basedir '/glmtools_misc/']);
addpath([basedir '/nlfuns/']);
addpath([basedir '/glmtools_spline/']);


%% load data! 

clear
% load 20180828_GLM_CouplingFilter_SameRGCtype 
CELL_TYPE = input('Cell type? (ON or OFF) ' ,'s')
fps = input('frame per second? (10 or 25) ')

%NUM_SET = input('Set number? ')
%CELL_TYPE = 'ON'
%NUM_SET = 1
% NUM_SET = 2
% NUM_SET = 3
% NUM_SET = 4
% NUM_SET = 5


loadDataNet 


%% First, run STA independently for all the channels 
%W = 1; % Window in sec
%W = 0.8; % Window in sec
%W = 0.5; % Window in sec
W = 0.6; % Window in sec


clear STAall
clear STAs  % STAs will be stored here

for i = 1:length(channelNames)
    
    channelName = channelNames{i}
    spikeTime=eval(channelName);
    
    [STA, gridT] = calcSTAprestim(stim, A1a, spikeTime, W, fps);
    %save STA
    STAall{i} = STA;

    
    %figure(i)
    clf
    STA_max_var = find_STA_with_max_var(STA, gridT, width, height)
    
    STAs(:,i) = STA_max_var; % store STA with maximum variance (center)
    
    % save figure 
    set(gcf, 'paperposition', [0 0 12 10])
    set(gcf, 'papersize', [12 10])

    saveas(gcf, sprintf('%scell_%dHz_%s_STA.pdf', CELL_TYPE, fps, channelName))
    
end

% plot STAs of multiple cells
clf
plot(gridT, STAs, 'linewidth',2)
ylabel('stim intensity')
xlabel('t')
title('STAs with the largest variance')
box off

legend(channelNames, 'location', 'NW','Interpreter', 'none'); legend boxoff


set(gcf, 'paperposition', [0 0 6 4])
set(gcf, 'papersize', [6 4])

saveas(gcf, sprintf('%scell_%dHz_STAs.pdf', CELL_TYPE, fps))


%% Now, let's fit GLM!

%% To do that, need to convert spike time to spike train

N = length(channelNames); % number of neurons 

binStim = size(stim,1);
assert (binStim==length(A1a))
spikeTrain = zeros(binStim,N);

for n = 1:N
    
    channelName = channelNames{n}
    spikeTime=eval(channelName);
    
    
    % let's check for each stim time bin
    for i = 2:binStim
        t0 = A1a(i-1);
        t1 = A1a(i);

        % find stim time that occured during t0 and t1
        idx = find(spikeTime>t0 & spikeTime<=t1);

        if ~isempty(idx)
            %disp('found')
            spikeTrain(i,n) = 1;
        end

    end
end


clf
imagesc(spikeTrain', [0 1]); colormap gray

% % Simple code below does NOT work due to time jitter!!!
% % grab spikeTime during stim & convert to bin idx
% idx = (spikeTime >= A1a(1)) & (spikeTime < A1a(end));
% binIdx = round(spikeTime(idx)*fps);  
% %idx = ceil(spikeTime*fps) % does it matter?
% spikeTrain(binIdx) = 1;
% spikeTrain = spikeTrain(:);
% 
% % cut spike train length to be the same as stim time length
% T = size(StimChosen,1);
% if length(spikeTrain) > T
%     spikeTrain = spikeTrain(1:T);
% end


%% Set up basis for coupling filters

% Make basis for self-coupling term
% ihbasprs.ncols = 2; % number of basis vectors
% ihbasprs.hpeaks = [.1, .2]; % peak of 1st and last basis vector
% ihbasprs.ncols = 3; % number of basis vectors
% ihbasprs.hpeaks = [.1, .2, .3]; % peak of 1st and last basis vector

switch fps
    case 10   % for 10 Hz datasets        
        % param for self-coupling filters
        ihbasprs.ncols = 1; % number of basis vectors
        ihbasprs.hpeaks = 0.1; % peak of 1st and last basis vector
        
        %ihbasprs.b = .001;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.b = .02;  % scaling (smaller -> more logarithmic scaling)
        %ihbasprs.b = 0.1;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.absref = []; % absolute refractory period basis vector (optional)

        % param for cross-coupling filters
        ihbasprs2.ncols = 1;  % number of basis vectors
        ihbasprs2.hpeaks = [0.01, .1]; % put peak at 0.1s and "effective" 1st peak at 0
        %ihbasprs2.b = .001;  % smaller -> more logarithmic scaling
        %ihbasprs2.b = .01;  % smaller -> more logarithmic scaling
        ihbasprs2.b = .1;  % smaller -> more logarithmic scaling
        ihbasprs2.absref = []; % no abs-refracotry period for this one

    case 25
        % param for self-coupling filters
        ihbasprs.ncols = 1; % number of basis vectors
        ihbasprs.hpeaks = 0.04; % peak of 1st and last basis vector
        
        %ihbasprs.b = .001;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.b = .02;  % scaling (smaller -> more logarithmic scaling)
        %ihbasprs.b = 0.1;  % scaling (smaller -> more logarithmic scaling)
        ihbasprs.absref = []; % absolute refractory period basis vector (optional)

        % param for cross-coupling filters
        ihbasprs2.ncols = 1;  % number of basis vectors
        ihbasprs2.hpeaks = [0.01, 0.04*2]; % put peak at 0.1s and "effective" 1st peak at 0
        %ihbasprs2.b = .001;  % smaller -> more logarithmic scaling
        %ihbasprs2.b = .01;  % smaller -> more logarithmic scaling
        ihbasprs2.b = .1;  % smaller -> more logarithmic scaling
        ihbasprs2.absref = []; % no abs-refracotry period for this one
        
end

% Make basis for self-coupling filter
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dtSpike);
nht = length(iht); % number of bins

% Make basis for cross-coupling filter
[iht2,ihbas2,ihbasis2] = makeBasis_PostSpike(ihbasprs2,dtSpike);
nht2 = length(iht2);

% pad to put them into the same time bins
if nht2>nht
    % padd ih1 with zeros
    iht = iht2; zz = zeros(nht2-nht,ihbasprs.ncols);
    ihbas = [ihbas;zz]; ihbasis = [ihbasis;zz]; nht=nht2;
elseif nht2<nht
    % padd ih1 with zeros
    iht2 = iht; zz = zeros(nht-nht2,ihbasprs2.ncols);
    ihbas2 = [ihbas2;zz]; ihbasis2 = [ihbasis2;zz]; nht2=nht;
end    


% Set an example filter shapes

%wself = [-10; 3.2; -1]; % weights for self-coupling term
wself = [-3; 0.5];
wself = wself(1:ihbasprs.ncols);
ihself = ihbasis*wself; % self-coupling filter
wcpl = 0.5; % weights for cross-coupling term
ihcpl = ihbasis2*wcpl; % cross-coupling filter


% plot coupling filter basis
%clf; 

c = get(0, 'DefaultAxesColorOrder');

XLIM = [0 0.7];
clf
subplot(311); 
plot(iht, ihbasis, 'Color', c(1,:), 'linewidth', 2); hold on
plot(ihbasprs.hpeaks, 1, '*k', 'linewidth', 1) %, 'MarkerSize', 7);
xlabel('time after spike (s)'); title('post-spike basis'); box off
set(gca,'xlim',XLIM);
subplot(312); 
plot(iht2, ihbasis2, 'Color', c(2,:), 'linewidth', 2); hold on
plot(ihbasprs2.hpeaks, 1, '*k', 'linewidth', 1) %, 'MarkerSize', 7);
xlabel('time after spike (s)'); title('coupling basis'); box off
set(gca,'xlim',XLIM);
subplot(313); 
plot(iht, exp(ihself), iht, exp(ihcpl), iht, iht*0+1, 'k:', 'linewidth', 2);
legend('self-coupling', 'cross-coupling');
xlabel('time lag (s)');
ylabel('gain (sp/s)');
box off
set(gca,'xlim',XLIM);

set(gcf, 'paperposition', [0 0 6 9])
set(gcf, 'papersize', [6 9])
saveas(gcf, sprintf('coupling-filters-basis_%dHz.pdf',fps))
saveas(gcf, sprintf('coupling-filters-basis_%dHz.png',fps))


%% Fit GLM 
%  Initialize params for fitting --------------
% Set params for fitting, including bases 


%ggInit = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

nkt  = size(STA,1);  % number of bins for stimulus filter
nkbasis = 2
ggInit = makeFittingStruct_GLM(dtStim,dtSpike,nkt,nkbasis,STA)


% Initialize fields for coupling filters (using h bases computed above)
ggInit.ihbas = ihbas; % h self-coupling basis
ggInit.ihbas2 = ihbas2; % h coupling-filter basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
nhbasis2 = size(ihbas2,2); % number of basis vectors in h basis
ggInit.ihw = zeros(nhbasis,1); % init params for self-coupling filter
ggInit.ihw2 = zeros(nhbasis2,N-1); % init params for cross-coupling filter
ggInit.ih = [ggInit.ihbas*ggInit.ihw ggInit.ihbas2*ggInit.ihw2];
ggInit.iht = iht;
ggInit.dc = 0; % Initialize dc term to zero
ggInit.sps = spikeTrain(:,1); % spikes from 1 cell
ggInit.sps2 = spikeTrain(:,2:N); % spikes from all 
ggInit.couplednums = 2:N; % cell numbers of cells coupled to this one 


ggfit(1:N) = ggInit; % initialize fitting struct
kfit = zeros(nkt,64,N);
dcfit = zeros(N,1);
for jj = 1:N
    couplednums = setdiff(1:N,jj);  % cell numbers of cells coupled to this one (jj)

    % Set spike responses for this cell (jj) and coupled cells
    ggInit.sps = spikeTrain(:,jj);
    ggInit.sps2 = spikeTrain(:,couplednums);
    ggInit.couplednums = couplednums; % numbers of cells coupled to this one (for clarity)

    % Do ML fitting
    fprintf('==== Fitting filters to neuron %d ==== \n',  jj)
    opts = {'display', 'iter', 'maxiter', 1000};
    ggfit(jj) = MLfit_GLM(ggInit,stim,opts); % do ML fitting
    kfit(:,:,jj) = ggfit(jj).k;
    dcfit(jj) = ggfit(jj).dc;    
end


%% plot linear filter of the last cell (for check)
clf;colormap gray
for n=1:N
    %%
    subplot(2,N,n);imagesc(STAall{n}');  ylabel('pixel'); xlabel('time bin');
    title(sprintf('STA of %s', channelNames{n}),'Interpreter', 'none')
    subplot(2,N,n+N);imagesc(ggfit(n).k'); ylabel('pixel'); xlabel('time bin')
    title(sprintf('Linear filter by GLM of %s', channelNames{n}),'Interpreter', 'none')
end
set(gcf, 'paperposition', [0 0 6 8])
set(gcf, 'papersize', [6 8])
saveas(gcf, sprintf('%scell_%dHz_STA_vs_GLM.pdf', CELL_TYPE,fps))


%% plot coupling filters
clf reset; 
colors = get(gca,'colororder');
ncolrs = min(size(colors,1),N); % number of traces to plot for each cell 

% find max of y range
ymax = 0;
for n=1:N
    ymax = max([exp(ggfit(n).ih(:)); ymax]); 
end
ymax = min(ymax,10);



switch N
    case 7 
        num_row=3;num_col=3;
    otherwise
        num_row=1;num_col=N;
end
    
for jj = 1:N
    subplot(num_row,num_col,jj); 
    plot(ggfit(jj).iht, exp(ggfit(jj).ih(:,1)), '--k', 'linewidth', 2); hold on
    plot(ggfit(jj).iht, exp(ggfit(jj).ih(:,2:ncolrs)), '-', 'linewidth', 2); % estim
    title(sprintf(' cell %s',channelNames{jj}(4:end) )); axis tight; set(gca,'ylim',[0,ymax]);
    box off
    
    if jj==1 || jj==4 || jj ==7
        ylabel('gain (sp/s)'); 
    end
    xlabel('time after spike (s)');
end


set(gcf, 'paperposition', [0 0 3*num_col 3*num_row])
set(gcf, 'papersize', [3*num_col 3*num_row])
saveas(gcf, sprintf('%scell_%dHz_coupling_filters.pdf', CELL_TYPE,fps))

