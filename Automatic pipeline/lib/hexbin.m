function [inbin, bincenters, vx, vy] = hexbin(pos, radius)
% [inbin, bincenters, vx, vy] = hexbin(pos, radius)
% hexagonal binning
% Xiaoyan, 2017

%% hexagon center and vertex coordinates
% distance between centers
dy = 3/2 * radius;
dx = sqrt(3) * radius;

% number of bins
nx = ceil(range(pos(:,1))/dx + .5);
ny = ceil(range(pos(:,2))/dy + .5);

% check if radius is too big
assert(nx>1, 'Hexbin radius is too big relative to spot x coordinates');
assert(ny>1, 'Hexbin radius is too big relative to spot y coordinates');

nx = nx + 1;
ny = ny + 1;

% X and Y center positions
cx = (min(pos(:,1))-dx/2) : dx : (min(pos(:,1))-dx/2+dx*nx);
cy = (min(pos(:,2))-dy/2) : dy : (min(pos(:,2))-dx/2+dy*ny);

% hexagon centers
[cx, cy] = meshgrid(cx, cy);

% shift x centers in every second horizontal lines
cx(1:2:ny,:) = cx(1:2:ny,:) + dx/2;

% vectorize and sort centers
centers = sortrows([cx(:), cy(:)]);
cx = centers(:,1);
cy = centers(:,2);

% vertex coordinates
[vy, vx] = dot2poly(cy, cx, radius, 6);


%% bin counts
inbin = zeros(size(pos,1),1);
counted = false(size(pos,1),1);

% count spots within a hexagon
for i = 1:length(cx)
    in = inpolygon(pos(:,1), pos(:,2), vx(:,i), vy(:,i));
    in = in & ~counted;
    counted = counted | in;
    inbin(in) = i;
end

% remove empty bins
bincenters = centers(unique(inbin),:);
vx = vx(:,unique(inbin));
vy = vy(:,unique(inbin));

end
