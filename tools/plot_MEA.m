function plot_MEA(offset_x, offset_y, LINE_TYPE)


%% default params
if nargin<1
    offset_x = 7;
end
if nargin<2
    offset_y = 4;
end 
if nargin<3
    LINE_TYPE='k';       %  default color is black
end
%hold on

grid_x = offset_x + 3*(0:7); % 전극간 거리는 좌표 ‘3’ 차이
grid_y = offset_y + 3*(0:7); % 전극간 거리는 좌표 ‘3’ 차이
[X,Y]=meshgrid(grid_x,grid_y);

%% mark circles
scatter(X(:), Y(:), LINE_TYPE)


%% mark channel numbers
hold on
% 
channel_no = 11:10:81;

text(grid_x-0.4, max(Y)+1, sprintfc('%d',channel_no))


axis xy

return 


%% To use
clf
plot_MEA