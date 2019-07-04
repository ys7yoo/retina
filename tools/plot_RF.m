function plot_RF(RF, FLIP_XY)

if nargin<2
    FLIP_XY = false;
end

if isempty(RF)
    return
end

if isfield(RF, 'type')
    switch RF.type
        case 'ON'
            LINE_STYLE = 'r-';
            TEXT_COLOR = [1 0 0];
        case 'OFF'
            LINE_STYLE = 'b-';
            TEXT_COLOR = [0 0 1];
    end
else
    % default line style
    LINE_STYLE = 'g-';
    TEXT_COLOR = [0 1 0];    
end


plot_ellipse(RF.mean, RF.cov, LINE_STYLE, FLIP_XY);

if isfield(RF, 'channel_name')
    if FLIP_XY
        tt=text(RF.mean(2), RF.mean(1), RF.channel_name(4:end), 'HorizontalAlignment','center');
    else
        tt=text(RF.mean(1), RF.mean(2), RF.channel_name(4:end), 'HorizontalAlignment','center');
    end
    tt.Color = TEXT_COLOR;
end    
        