
%% load data
clear
% set 1
% load('20180404/StimInfo_Num1_13pix_100um_10Hz.mat')
% fps = 10
% set 2 
load('20180404/StimInfo_Num3_13pix_100um_20Hz.mat')
fps = 20
height=13
width=13

T = length(StimInfo)

stim = cell2mat(StimInfo);
stim = reshape(stim, [T, height, width]);
% rearrange order
stim = permute(stim, [2 3 1]);
size(stim)

% plot some sample stimulus
clf
subplot(131)
imshow(stim(:,:,1))
subplot(132)
imshow(stim(:,:,2))
subplot(133)
imshow(stim(:,:,3))


%% Let's analyze one cell from Dataset 1
% ON cell (ch47b,58a,58b,66c,66d,76a,77b,85c,86c)
% OFF cell (ch16b,17a,17b,28a,37a,46a,47a,48b,55b,65b,66a,66b,67a,68b,77a,85a,86a)

% % ch47b
% spikeTime = [17.74784
% 27.6074
% 31.75488
% 32.69296
% 35.07048
% 35.1534
% 35.18752
% 35.2962
% 36.89488
% 37.16388
% 38.55844
% 44.1116
% 44.16504
% 44.83928
% 45.5786
% 47.20044
% 47.21612
% 54.67052
% 55.67812
% 58.1632
% 58.79664
% 63.03024
% 67.23396
% 69.1698
% 69.19332
% 69.24792
% 69.52464
% 73.99512
% 77.06608
% 77.07144
% 77.2986
% 84.42076
% 86.43788
% 90.57816
% 90.72556
% 93.1456
% 93.16776
% 97.6384
% 97.884
% 101.39348
% 103.70556
% 103.74248
% 104.39224
% 106.1778
% 115.28372
% 116.02552
% 119.0414
% 120.53984
% 121.40572
% 124.94188
% 125.38756
% 135.94512
% 139.95776
% 141.62404
% 142.23864
% 146.30236
% 149.13
% 150.89956
% 152.0442
% 153.67384
% 154.63004
% 155.43288
% 159.31732
% 160.52828
% 161.65004
% 162.60512
% 164.65196
% 169.6246
% 176.123
% 176.63796
% 176.65904
% 177.2244
% 181.62156
% 183.8906
% 190.42004
% 192.96848
% 200.39
% 201.89572
% 204.37984
% 209.57244];


% from dataset 3 (20Hz)

spikeTime = [4.26488
4.6826
5.62696
6.42896
7.53348
9.54728
10.26096
13.58896
13.6154
14.48352
15.27024
17.67932
22.31172
22.69444
23.1886
25.06608
29.05508
29.28068
31.78316
32.34056
35.43216
35.74896
36.28388
38.67268
39.14584
39.28468
46.3644
51.20904
51.41004
51.68008
52.47864
52.80564
52.8628
55.88748
56.39416
57.38436
59.56476
64.03596
64.26368
64.5624
66.53544
67.01032
67.01628
69.57184
69.77032
70.89756
74.0742
74.283
77.3118
77.5232
77.68104
77.95944
81.70612
84.0708
85.54728
85.84272
86.5054
86.5228
86.73476
88.2874
88.50636
93.341
96.70748
97.40084
99.98528
100.21492
100.4832
100.59584
101.03956
101.87076
102.5774
104.14176
105.81148
106.01324
106.75628
109.00948
112.75604
114.6856
115.51688
115.74176
118.33268
119.8562
120.11512
122.23908
122.9986
128.44388
128.62888
128.67676
129.16556
129.19668
129.92224
130.04876
130.0784
130.74332
131.16
131.18024
133.47456
136.29676
137.46884
139.94992
140.33504
140.76332
142.56084
143.41644
147.20888
148.58084
149.33064
149.59088
149.86684
150.1132
150.39096
150.64128
153.20584
154.73444
160.68192
161.03604
161.47692
161.65272
162.51192
162.97056
163.20468
165.1238
165.2734
165.41676
166.88988
174.45068
179.9132
180.57436
181.24336
181.88148
182.38868
183.61448
184.0642
188.12748
188.73728
190.28852
191.24444
197.9116
199.16172
200.67204
200.90252
202.03424
202.46616
202.8046
204.96988
205.24332
207.96288
207.98144
208.28568
208.62208
209.19904];



% convert to index
spikeTimeIdx = round(spikeTime*fps)

if (max(spikeTimeIdx) > T)
    warning('Maximum spike time exceeds recoding time. Check!!!')
    
    % [QUICK FIX]
    % remove spikeTimeIdx out of range (recording time shift???) [TO CHECK]
    spikeTimeIdx = spikeTimeIdx(spikeTimeIdx<=T);
end
%%  calc STA

% Step 1. find stimulus pattern that triggers a spike
W = 1*fps;    % window size 
stimTriggeredSpike = [];
cnt = 1;

% remove spikeTimeIdx less than W
spikeTimeIdx = spikeTimeIdx(spikeTimeIdx>W);

numSpike = length(spikeTimeIdx);

stimTriggeredSpike = zeros (height,width,W,numSpike);
for i=1:length(spikeTimeIdx)
    idx= spikeTimeIdx(i);
        stimTriggeredSpike(:,:,:,i) = stim(:,:,idx-W+1:idx);
end

STA = mean(stimTriggeredSpike,4); % 3-dim array height x width x W

%% plot STA
STAstack = reshape(STA,[height*width, W]);
clf
subplot(121)
imshow(STAstack)
xlabel('t')
ylabel('pixel')
title('STA')

% check variance to see if there is any change 
STAvar = var(STAstack,[],2);
subplot(222)
imshow(reshape(STAvar,height,width),[])
xlabel('x')
ylabel('y')
title('variance across time')


% plot STA for the pixel with the largest variance
[mm, maxIdx] = max(STAvar)
subplot(224)
plot(STAstack(maxIdx,:))
xlabel('t')
ylabel('STA')
title('STA for the pixel with the largest variance')


set(gcf, 'paperposition', [0 0 9 8])
set(gcf, 'papersize', [9 8])
saveas(gcf, 'STA.pdf')

                                 