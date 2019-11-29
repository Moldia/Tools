function density = gene_kde(name, pos, query, bandwidth, imfile, scale)
% density = gene_kde(name, pos, query, bandwidth, imfile, scale)
% density estimate of a gene
% output is at image scale/5 
% Xiaoyan, 2018

% image size
imsize = imfinfo(imfile);
imsize = [imsize.Height, imsize.Width];

% scale
if nargin <= 5
    scale = 1;
end

pos = correctcoord(pos, .2);
[uNames, ~, iName] = unique(name);

iNamePlot = find(strcmp(uNames, query));

if ~isempty(iNamePlot)
    posPlot = pos(iName == iNamePlot,:);
    try
        [~, density] = kde2d_modified...
            (posPlot, 2^10, [0 0],...
            floor([imsize(2)/scale/5, imsize(1)/scale/5]),...
            floor([bandwidth/5, bandwidth/5]));
    catch
        density = zeros(2^10);
        warning('Too few transcripts');
    end
else
    warning('No specified transcript detected in the input file');
    density = zeros(2^10);
end
density = imresize(density, [imsize(1)/5 imsize(2)/5]);
end
