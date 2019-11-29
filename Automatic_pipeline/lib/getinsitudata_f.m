function [name,pos,quality,cell] = getinsitudata_f(filename,varargin)

% Return the names and positions from the input file
% Last update, Xiaoyan 2016-4-17

cell_column = 0;
switch length(varargin)
    case 0
        name_column = 2;
        pos_column = 1;
        score_column = 5;
    case 1
        name_column = varargin{1};
        pos_column = 1;
        score_column = 5;
    case 2
        name_column = varargin{1};
        pos_column = varargin{2};
        score_column = 5;
    case 3
        name_column = varargin{1};
        pos_column = varargin{2};
        score_column = varargin{3};
    case 4
        name_column = varargin{1};
        pos_column = varargin{2};
        score_column = varargin{3};
        cell_column = varargin{4};
end

data = importdata(filename,',',1);
name = data.textdata(2:end,name_column);
pos = data.data(:,pos_column:pos_column+1);
if score_column
    quality = data.data(:,score_column);
end
if cell_column
    cell = data.data(:,cell_column);
end
end