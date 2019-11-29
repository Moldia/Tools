function [uNames, cName] = plotall(name, pos, imfile, scale, second_imfile)
% [uNames, cName] = plotall(name, pos, imfile, scale, second_imfile)
% simply plot all reads
% Xiaoyan, 2018

% prepare background
% figure;

if ischar(imfile)
    try
        imfile = imread(imfile);
        imshow(imfile, []);
    catch
        plotonblank;
    end
else
    imshow(imfile, []);
end

if nargin <= 3
    scale = 1;
elseif nargin > 4
    try
        addimlayer(second_imfile, .4);
    end
end

% unique genes
[uNames, ~, iName] = unique(name);
cName = hist(iName, 1:length(uNames));

pos = correctcoord(pos, scale);
sym = repmat(symlist, ceil(length(uNames)/length(symlist)), 1);
hold on
for i = 1:length(uNames)
    plot(pos(iName==i,1), pos(iName==i,2), sym{i});
end

legend(uNames, 'color', [.6 .6 .6], 'location', 'NorthEastOutside');

end
