function [codebook, plProbes, barcodeLetter, swProbes] = ...
    get_codebook(query_oligonum, db_oligonum, db_genename, db_barcode,...
    seqorder, base_channel_order, DO_base)
% [codebook, plProbes, barcodeLetter, swProbes] = ...
%     get_codebook(query_oligonum, db_oligonum, db_genename, db_barcode,...
%     seqorder, base_channel_order, DO_base)
% prepare codebook
% Xiaoyan, 2017

% find probes
hits = cellfun(@(v) strcmp(v, db_oligonum), query_oligonum, 'uni', 0);
hits = cellfun(@find, hits);

% separate sandwich probes
swProbes = cellfun(@(v) strcmp(v(1:2), 'SW'), db_barcode(hits));
plProbes = hits(~swProbes);
swProbes = hits(swProbes);

% convert DO to bases
if nargin > 6
    barcodeLetter = cellfun(@(v) [DO_base{str2double(v(1))}, v(2:end)],...
        db_barcode(plProbes), 'uni', 0);
else
    barcodeLetter = db_barcode(plProbes);
end

% base to number
[~, barcodenum] = barcode2num(barcodeLetter);

% sequencing order
if nargin < 5
    seqorder_numeric = 1:length(barcodeLetter{1});
else
    seqorder_numeric = [];
    for i = 1:length(seqorder)
        seqorder_numeric(i) = str2double(seqorder(i));
    end
end

% adjust based on sequencing order
barcodeLetter = cellfun(@(v) reorder_letters(v, seqorder_numeric),...
    barcodeLetter, 'uni', 0);
barcodeDigit = barcodenum(:,seqorder_numeric);

% convert to channel number
barcodeChannel = base_channel_order(barcodeDigit);
barcodeChannel = barcodeChannel +...
    repmat(0:6:6*(length(seqorder_numeric)-1), length(plProbes), 1); % in full dimension

% add SW probes to the end
SWchannels.Sst = 1;
SWchannels.Npy = 3;
if ~isempty(swProbes)
    barcodeChannel = [barcodeChannel, zeros(size(barcodeChannel,1),1)];
end
for i = 1:length(swProbes)
    barcodeChannel = [barcodeChannel;...
        [zeros(1,size(barcodeChannel,2)-1),...
        6*length(seqorder_numeric) + base_channel_order(SWchannels.(db_genename{swProbes(i)}))]];
end

% binary code matrix
barcodeChannel(barcodeChannel==0) = 6*(length(seqorder_numeric)+nnz(~isempty(swProbes)))+1;
barcodeDigit = repmat((1:length(query_oligonum))', 1, size(barcodeChannel,2));
codebook = accumarray([barcodeDigit(:),barcodeChannel(:)],...
    1, [length(query_oligonum), max(barcodeChannel(:))]);
if ismember(0, barcodeChannel)
    codebook(:,end) = [];   % remove last base N
end

% general stain
codebook(1:end-length(swProbes), 2:6:6*length(seqorder_numeric)) = 1;

end


function barcode_new = reorder_letters(barcode, order)
barcode_new = char('');
for i = 1:length(order)
    barcode_new(i) = barcode(order(i));
end
end
