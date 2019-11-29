function imLabel = relabelimg(imLabel)
% imLabel = relabelimg(imLabel)
% relabel label images (start with label=1) outputed from Cellprofiler
% Xiaoyan, 2017

tilesize = size(imLabel,1);
idxLabel = 0;
imLabel = imLabel(:);
uniLabel = unique(imLabel);

if length(uniLabel) == 1
    imLabel(:) = 0;
else
    % convert 16-bit numbers to object numbers
    for i = 1:length(uniLabel)
        imLabel(imLabel==uniLabel(i)) = idxLabel;
        idxLabel = idxLabel+1;
    end
end

imLabel = reshape(imLabel,tilesize,[]);
end
