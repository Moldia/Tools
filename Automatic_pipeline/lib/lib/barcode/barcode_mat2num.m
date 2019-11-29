function barcode_num = barcode_mat2num(barcode_matrix)
% barcode_num = bardcode_mat2num(barcode_matrix)
% convert barcode matrix to number list
% Xiaoyan, 2017


nBase = size(barcode_matrix, 2);
barcode_num = zeros(size(barcode_matrix, 1), 1);

for b = 1:nBase
    barcode_num = barcode_num + barcode_matrix(:,b)*10^(nBase-b);
end

