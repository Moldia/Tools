function [poly_x, poly_y] = dot2poly(x, y, radius, nvertices)
% [poly_x, poly_y] = dot2poly(x, y, radius, corner)
% Convert dot to polygon circles
% Xiaoyan 2014-11-18



[circx, circy] = circlepolygon(radius, nvertices);
poly_x = (repmat(x,1,nvertices+1) + repmat(circx,length(x),1))';
poly_y = (repmat(y,1,nvertices+1) + repmat(circy,length(y),1))';


    function [circx, circy] = circlepolygon(radius, nvertices)
        theta = linspace(0, 2*pi, nvertices+1);
        circx = radius*cos(theta);
        circy = radius*sin(theta);
    end

end
