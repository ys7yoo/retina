function plot_RF(RF, FLIP_XY)

if nargin<2
    FLIP_XY = false;
end

if ~isfield(RF, 'type')
    return
end

switch RF.type
    case 'ON'
        LINE_STYLE = 'r-';
        TEXT_COLOR = [1 0 0];
    case 'OFF'
        LINE_STYLE = 'b-';
        TEXT_COLOR = [0 0 1];        
end


plot_ellipse(RF.mean, RF.cov, LINE_STYLE, FLIP_XY);
if FLIP_XY
    tt=text(RF.mean(2), RF.mean(1), RF.channel_name(4:end), 'HorizontalAlignment','center');
else
    tt=text(RF.mean(1), RF.mean(2), RF.channel_name(4:end), 'HorizontalAlignment','center');
end
tt.Color = TEXT_COLOR;
        