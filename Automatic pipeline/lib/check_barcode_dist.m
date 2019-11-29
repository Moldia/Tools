function check_barcode_dist(barcodefile, nHybs, isDO, varargin)
% check_barcode_dist(barcodefile, nHybs, isDO, varargin)
% input: barcode file(csv) - column 1: name, column 2: barcode
% output: distance matrix figure and probe pairs with too close distance
% Xiaoyan, 2017


barcodelist = importdata(barcodefile);
barcodelist = cellfun(@(v) strsplit(v, ','), barcodelist, 'uni', 0);
barcodelist = reshape([barcodelist{:}], 2, [])';

% take only normal plp
barcodes = barcodelist(:,2);
lenbarcode = cellfun(@length, barcodes);
normplp = lenbarcode == nHybs+isDO;
name = barcodelist(normplp,1);

if isDO    % 1st letter in barcode is numeric, repsresenting detection oligo number
    barcodes = barcodes(normplp);
    DO = cell2mat(cellfun(@(v) str2num(v(1)), barcodes, 'uni', 0));
    barcodes = cellfun(@(v) v(2:end), barcodes, 'uni', 0);
else
    barcodes = barcodes(normplp);
end

% change to numbers
try
    [~, nummatrix] = barcode2num(barcodes, varargin{1});
catch
    [~, nummatrix] = barcode2num(barcodes);
end
if isDO
    nummatrix = [nummatrix, DO];
end

% barcode distance
Dist = zeros(size(nummatrix,1), 'uint8');
for i = 1:(nHybs + isDO)
    Dist = Dist + uint8(...
        repmat(nummatrix(:,i), 1, size(nummatrix,1)) ~=...
        repmat(nummatrix(:,i)', size(nummatrix,1), 1));
end
Dist(Dist>2) = 2;

% figure
figure,imshow(Dist,[]);
colormap hot
axis on

% close pairs
[m, n] = find(Dist==1);
lowerhalf = find(m>n);
for i = 1:length(lowerhalf)
    fprintf('(%s, %s)\n', name{m(lowerhalf(i))}, name{n(lowerhalf(i))});
end

end
