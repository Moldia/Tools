function [gridorder, imfiles] = tilepos2grid(tilepos, imprefix, imsuffix, ndigits)
% convert metadata tile position to grid
% Xiaoyan, 2017


% find stage movement delta pixel
tilepos = round(tilepos, -2);
tilesize = unique(tilepos(:));
tilesize = abs(repmat(tilesize, 1, length(tilesize)) -...
    repmat(tilesize', length(tilesize), 1));
[ndist, tilesize] = hist(tilepos(:), unique(tilepos(:)));
% [~, idx] = sort(ndist, 'descend');
tilesize = tilesize(2);

% convert to grid
gridpos = floor(tilepos/tilesize) + 1;
ntiles = max(gridpos, [], 1);

gridorder = zeros(ntiles);
gridorder((gridpos(:,2)-1)*ntiles(1) + gridpos(:,1)) = 1:length(tilepos);

% almost all microscope stages move along x axis first
gridorder = gridorder'; 

imfiles = cell(size(gridorder));
if nargin > 1
	try
		imfiles = catstrnum(imprefix, gridorder, ndigits);
    catch
		imfiles = catstrnum(imprefix, gridorder);
	end
    imfiles = strcat(imfiles, imsuffix);
    imfiles(gridorder==0) = {'empty.tif'};
end


end
