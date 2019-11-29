function I = resizestitch(ntiles, tilesize, resizef, immatrix)
% resize and stitch tile images
% Xiaoyan, 2018


if length(size(immatrix{1}))==3
    I = zeros(ntiles(2)*tilesize*resizef, ntiles(1)*tilesize*resizef, 3, 'uint8');
else
    I = zeros(ntiles(2)*tilesize*resizef, ntiles(1)*tilesize, 'uint32');
end

for i = 1:ntiles(2)
    for j = 1:ntiles(1)
        temp = imresize(immatrix{i,j}(:,:,:),resizef);
        I(tilesize*resizef*(i-1)+1:tilesize*resizef*i,...
            tilesize*resizef*(j-1)+1:tilesize*resizef*j,:) = temp;
    end
end

end
