%% set path to GLM tools

% This is where the GLMspiketools is installed
% You can get the package from https://github.com/ys7yoo/GLMspiketools
basedir = [getenv('HOME') '/src/GLMspiketools']

% Add a bunch sub-directories (with absoluate path names)
addpath([basedir '/glmtools_fitting/']);
addpath([basedir '/glmtools_misc/']);
addpath([basedir '/nlfuns/']);
addpath([basedir '/glmtools_spline/']);


%% First run STA
runSTAprestim

%% choose pixel to use 
clf
subplot(211)
imagesc(STA)
xlabel('pixel index')
ylabel('time bin')
axis xy

subplot(212)
plot(STAvar)
xlabel('pixel index')
ylabel('variance over time')



%% choose first 8 pixels (for simple test)

idxChannel = 1:8
%idxChannel = 1:5
STAchosen = STA(:,idxChannel);
StimChosen = stim(:,idxChannel);

clf
subplot(211)
imagesc(StimChosen')
title('Stim for chosen pixels')
ylabel('chosen pixel index')
xlabel('time bin (in stim dT)')

subplot(212)
imagesc(STAchosen')
title('STA for chosen pixels')
ylabel('chosen pixel index')
xlabel('time bin (in spike dT)')
axis xy

%% Now, let's fit GLM!

%% To do that, need to convert spike time to spike train

binStim = size(StimChosen,1);
assert (binStim==length(A1a))

spikeTrain = zeros(binStim,1);
% let's check for each stim time bin
for i = 2:binStim
    t0 = A1a(i-1);
    t1 = A1a(i);
    
    % find stim time that occured during t0 and t1
    idx = find(spikeTime>t0 & spikeTime<=t1);
    
    if ~isempty(idx)
        %disp('found')
        spikeTrain(i) = 1;
    end
    
end

clf
plot(spikeTrain)

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



%% Fit GLM 
%  Initialize params for fitting --------------
% Set params for fitting, including bases 
dtStim = 1 / fps;
dtSpike = 1 / fps;
nkt  = size(STAchosen,1);  % number of bins for stimulus filter
nkbasis = 4;  % number of basis vectors for representing k
nhbasis = 2;  % number of basis vectors for representing h
hpeakFinal = .1;   % time of peak of last basis vector for h
gg0 = makeFittingStruct_GLM(dtStim,dtSpike,nkt,nkbasis,STAchosen)
%gg0 = makeFittingStruct_GLM(dtStim,dtSpike,nkt,nkbasis,STAchosen,nhbasis,hpeakFinal);


gg0.sps = spikeTrain;  % Insert binned spike train into fitting struct
gg0.mask = [];  % Not currently supported!; % insert mask (optional)
gg0.ihw = randn(size(gg0.ihw))*1; % initialize spike-history weights randomly

% Compute conditional intensity at initial parameters 
[negloglival0,rr] = neglogli_GLM(gg0,StimChosen);
fprintf('Initial negative log-likelihood: %.5f\n', negloglival0);

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg1, negloglival1a, H1, Xstruct1] = MLfit_GLM(gg0,StimChosen,opts); % do ML (requires optimization toolbox)



%% Fit GLM ("bilinear stim filter version") via max likelihood

%  Initialize params for fitting --------------
k_rank = 1; % Number of column/row vector pairs to use
gg0b = makeFittingStruct_GLMbi(k_rank,dtStim,dtSpike,nkt,nkbasis,STAchosen);
%gg0b = makeFittingStruct_GLMbi(k_rank,dtStim,dtSpike,nkt,nkbasis,STAchosen,nhbasis,hpeakFinal);
gg0b.sps = spikeTrain;
gg0b.mask = []; 
logli0b = neglogli_GLM(gg0b,StimChosen); % Compute logli of initial params
fprintf('Initial value of negative log-li (GLMbi): %.3f\n', logli0b);

% Do ML estimation of model params
opts = {'display', 'iter'};
[gg2, negloglival2] = MLfit_GLMbi(gg0b,StimChosen,opts); % do ML (requires optimization toolbox)


%%  plot results
clf
subplot(321)
imagesc(STAchosen')
title('STA')
ylabel('chosen pixel index')
xlabel('time bin')
axis xy

subplot(322)
plot(STAchosen)
xlabel('time bin')
axis tight
box off

subplot(323)
imagesc(gg1.k'); title('Estimate of K by full GLM'); ylabel('chosen pixel index'); xlabel('time bin');
axis xy

subplot(324)
plot(gg1.k)
xlabel('time bin')
axis tight
box off

subplot(325)
imagesc(gg2.k'); title('Estimate of K by bi-linear GLM'); ylabel('chosen pixel index'); xlabel('time bin');
axis xy

subplot(326)
plot(gg2.k)
xlabel('time bin')
axis tight
box off

% save figure
set(gcf, 'paperposition', [0 0 10 8])
set(gcf, 'papersize', [10 8])
saveas(gcf, sprintf('Num%d_GLM_%s.pdf',NUM_EXP,channelName))


return

%% Let's further understand the effect of each component


% gg1, negloglival1a, H1, Xstruct1

prs1 = [gg1.kt(:); gg1.dc; gg1.ihw; gg1.ihw2];
Loss_GLM_logli_exp(prs1,Xstruct1)



