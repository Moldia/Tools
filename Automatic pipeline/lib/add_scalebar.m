function add_scalebar(pixelscaling, unit, barpos, col, textpos, fontsize, orientation)
% add_scalebar(pixelscaling, unit, barpos, col, textpos, fontsize, orientation)
% add scalebar to the current image
% all inputs are optional
% Xiaoyan, 2018


if nargin < 1
    pixelscaling = 0.33;
end

if nargin <= 1
    unit = 1000;
end

xrange = get(gca, 'xlim');
yrange = get(gca, 'ylim');

if nargin <= 2 || isempty(barpos)
    barpos(1) = (xrange(2) - unit/pixelscaling)*.95;
    barpos(2) = yrange(2)*.97;
end

if nargin <= 3
    col = 'w';
end

if nargin <= 4
    textpos = {'bottom', 'center'};
end

if nargin <= 5
    fontsize = 10;
end

if nargin <= 6
    orientation = 'horizontal';
end


hold on;

if strcmp(orientation, 'horizontal')
    plot([barpos(1), barpos(1)+unit/pixelscaling], [barpos(2) barpos(2)],...
        col, 'linewidth', 2);
else
    plot([barpos(1), barpos(1)], [barpos(2) barpos(2)-unit/pixelscaling],...
    col, 'linewidth', 2);
end

if ~isempty(textpos)
    if strcmp(textpos, 'default')
        textpos = {'bottom', 'center'};
    end
        
    if strcmp(orientation, 'horizontal')
        text(barpos(1)+unit/pixelscaling/2, barpos(2),...
            [num2str(unit) ' \mum'],...
            'verticalalignment', textpos{1}, 'horizontalalignment', textpos{2},...
            'color', col, 'fontsize', fontsize);
    else
        text(barpos(1), barpos(2)-unit/pixelscaling/2,...
            [num2str(unit) ' \mum'],...
            'verticalalignment', textpos{1}, 'horizontalalignment', textpos{2},...
            'color', col, 'fontsize', fontsize,...
            'Rotation', 90);
    end
end

end

