function [setNumber, setLetter] = barcode_generator(nHybs, nDiffBases,...
    OccupiedBarcodes, codedist, letters, maxstretch)
% [setNumber, setLetter] = barcode_generator(nHybs, nDiffBases,...
%     OccupiedBarcodes, codedist, letters)
% generate all possible barcodes that are certain distance away from each other
% OccupiedBarcodes (optional): a list of already used barcodes (letter or number)
% codedist (optional): barcode distance, default 2
% letters (optional): letter list, default {'A' 'C' 'G' 'T'}
% Xiaoyan, 2017


% all possible barcodes
[Y{nHybs:-1:1}] = ndgrid(1:nDiffBases);
nummatrix = reshape(cat(nHybs+1, Y{:}), [], nHybs);
numlist = zeros(size(nummatrix,1), 1);
for i = 1:nHybs
    numlist = numlist + nummatrix(:,i)*10^(nHybs-i);
end
clear Y

% remove barcodes that have same base stretching more than maxstretch
if nargin>5 && maxstretch<nHybs
    stretch = ~logical(nummatrix(:,2:end) - nummatrix(:,1:end-1));
    for i = 1:min(maxstretch-1, nHybs)
        stretch = stretch(:,2:end) & stretch(:,1:end-1);
    end
    stretch = logical(sum(stretch,2));
    numlist = numlist(~stretch);
    nummatrix = nummatrix(~stretch, :);
end

% homomer
homomatrix = repmat((1:nDiffBases)', 1, nHybs);
homomer = zeros(nDiffBases,1);
for i = 1:nHybs
    homomer = homomer + homomatrix(:,i)*10^(nHybs-i);
end
homomer = num2cell(homomer);
homomer = cellfun(@(v) v==numlist, homomer, 'uni', 0);
homomer = logical(sum([homomer{:}], 2));

% already used barcodes
if nargin < 5
    letters = {'A' 'C' 'G' 'T'};
end
if nargin > 2
    if isnumeric(OccupiedBarcodes)
        occupied = cellfun(@(v) v==numlist, num2cell(OccupiedBarcodes), 'uni', 0);     % varargin{1} is numeric
    else
        occupied = cellfun(@(v) v==numlist,...
            num2cell(barcode2num(OccupiedBarcodes, letters)), 'uni', 0);
    end
end

% convert into logical
if exist('occupied', 'var') && ~isempty(occupied)
    occupied = logical(sum([occupied{:}], 2));
else
    occupied = false(length(numlist),1);
end

% desired barcode distance
if nargin <= 3
    codedist = 2;
end

% barcode pairwise distance
disp('calculating pairwise distance..');
cellnummatrix = num2cell(nummatrix,2); 
iNN = cellfun(@(v) find(calculate_base_dist(v, nummatrix)<codedist), cellnummatrix,...
    'uni', 0);
% iNN = {};
% parfor i = 1:length(numlist)
%     iNN{i} = find(calculate_base_dist(nummatrix(i,:), nummatrix)<codedist);
% end

% neighbors of occupied
occupiedNeighbors = false(length(numlist),1);
occupiedNeighbors(cat(1, iNN{occupied})) = true;

% exhaustive search to find all possible barcode SETS
setNumber = [];
if nargout > 1
    setLetter = [];
end

disp('finding usable barcode sets..');
while nnz(~occupied)
    discard = false(length(numlist) ,1);
    discard = discard | occupied | homomer | occupiedNeighbors;
    
    for i = 1:length(numlist)
        if ~discard(i)
            % discard all neighbors of currently chosen one
            discard([iNN{i}]) = true;
            % but keep itself
            discard(i) = false;
        end
    end
    
    % the ones chosen (i.e. not discarded) in the current loop
    currentset = find(~discard);
    occupied = occupied | ~discard;
    if ~isempty(currentset)
        setNumber = [setNumber, {numlist(currentset)}];
    else
        break
    end

    try
        setLetter = [setLetter, {num2barcode(numlist(currentset), letters)}];      
    catch me
        if strfind(me.message, 'letters')
            error('Number of different bases exceeds default ACGT.');
        elseif strfind(me.identifier, 'Undefined')
            continue
        end
    end
end

    function d = calculate_base_dist(query, pool)
        d = zeros(length(numlist), 1, 'uint8');
        for b = 1:nHybs
            d = d + uint8(query(b)~=pool(:,b));
        end
    end
end

