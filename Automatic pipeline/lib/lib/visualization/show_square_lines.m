function show_square_lines(f, sizeX, sizeY, varargin)
% draw square lines, usu. bin or tile lines
% Xiaoyan, 2017


figure(f);
hold on;

if ~isempty(varargin)
    col = varargin{1};
else
    col = 'b';
end

if length(varargin)<=1
    offset = 0;
else
    offset = varargin{2};
end

ax = gca;
limX = ax.XLim;
limY = ax.YLim;

x = sizeX;
while x < limX(2)
    plot([x,x]+offset, limY+offset, col);
    x = x + sizeX;
end

y = sizeY;
while y < limY(2)
    plot(limX+offset, [y, y]+offset, col);
    y = y + sizeY;
end

end