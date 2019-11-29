function transformed = transformimage(original, ang, yup, xleft, imscale, name, refsize, saveimage)
% transformed = transformimage(original, ang, yup, xleft, imscale, name, refsize, saveimage)
% transform image and make it the same as as reference
% Xiaoyan, 2017

transformed = imrotate(original, ang);
if yup>=0
    transformed = transformed(yup*imscale+1:end,:);
else
    transformed = [zeros(-yup*imscale,size(transformed,2)); transformed];
end
if xleft>=0
    transformed = transformed(:,xleft*imscale+1:end);
else
    transformed = [zeros(size(transformed,1),-xleft*imscale), transformed];
end

% keep the transformed image the same size as specified (crop/ pad lower
% right corner)
if nargin >= 7
    size_delta = refsize-size(transformed);
    transformed = padimg(transformed, size_delta(2), size_delta(1));
end
    
% save if output filename is given
if nargin>5 && nargin<8 || saveimage
    imwrite(transformed, name, 'tiff', 'compression', 'none');
end

end
