function h = plot_MEA_param(ab, cd)

% plot MEA using fitted params

if nargin<1 
    ab = [3 0];  % default param for x
end

if nargin<2
    cd = [-3 0];
end

% if nargin<3
%     LINE_TYPE='k';       %  default color is black
% end
%%

grid_y = ab(1)*(1:8) + ab(2);   % 전극간 거리는 좌표 ‘3’ 차이
grid_x = cd(1)*(1:8) + cd(2);   % 전극간 거리는 좌표 ‘3’ 차이
[Y,X]=meshgrid(grid_y,grid_x);

%%
h = scatter(Y(:), X(:), [], [128, 128, 128]/255);


% %% mark channel numbers
% hold on
% % 
% channel_no = 11:10:81;
% 
% text(grid_y, max(X)+1.5, sprintfc('%d',channel_no), 'HorizontalAlignment','center')


axis xy
