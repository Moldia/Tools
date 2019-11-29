function addimlayer(imf, alpha)
% addimlayer(imf, alpha)
% overlay a transparent image layer
% Xiaoyan, 2017

hold on;
if ischar(imf)
    im = imread(imf);
else
    im = imf;
end
ih = imshow(im);
set(ih, 'alphadata', alpha);

end
