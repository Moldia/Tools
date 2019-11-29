function [codemat, codebook] = code2binarymat(taglist, base_channel_order, anchor_channel)
% [codemat, codebook] = code2binarymat(taglist, base_channel_order, anchor_channel)
% convert UNIQUE barcodes into binary matrix
% Xiaoyan, 2017


% correct for potential barcode duplicates
taglist = remove_barcode_duplicates(taglist);

% base to number
barcode = taglist(:,2);
[~, barcodeDigit] = barcode2num(barcode);

% channel number
barcodeChannel = base_channel_order(barcodeDigit);
barcodeChannel = barcodeChannel +...
    repmat(0:6:6*(length(barcode{1})-1), length(barcode), 1); % in full dimension


% binary code matrix
barcodeChannel(barcodeChannel==0) = 6*length(barcode{1})+1;
barcodeDigit = repmat(1:length(barcode)', 1, size(barcodeChannel,2));
codemat = accumarray([barcodeDigit(:),barcodeChannel(:)],...
    1, [length(barcode), max(barcodeChannel(:))]);

% general stain
codemat(:,anchor_channel:6:6*length(barcode{1})) = 1;

% ready to write format
header = {'gene', 'barcode'};
for i = 1:length(barcode{1})
    for j = 1:6
        header = [header, {strcat('base', num2str(i), '_channel', num2str(j))}];
    end
end
codebook = [header;...
    taglist(:,1), barcode, num2cell(codemat)];

end

