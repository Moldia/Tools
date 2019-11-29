function inROI = readsinroi(pos, imLabel)
% inROI = readsinroi(pos, imLabel)
% find reads in labeled mask(s)
% Xiaoyan, 2017

pos = round(pos);
posLinear = (pos(:,1)-1)*size(imLabel,1) + pos(:,2);
inROI = zeros(length(posLinear), 1);
for i = 1:max(double(imLabel(:)))
    inBW = imLabel(posLinear)==i & inROI==0;
    inROI(inBW) = i;
end
end
