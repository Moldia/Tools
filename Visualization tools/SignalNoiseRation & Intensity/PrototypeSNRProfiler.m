
image = '/Volumes/SBHvsSBL/imagesAndCSVFilesFor3ROI/Round4Cycle4/SBH_Round4Cycle4_SBH_ROI1/SBH_Round4Cycle4_SBH_ROI1_c2_ORG.tif';
  
I = imread(image);
maxdist_between_spots = 30;
filename = [strtok(image, '.'), '.csv']; % strok is used to select parts of strings
% we initiate a file with the file name that comes from image and we then
% append the proper .csv extension

%threshold = prctile(I(:), 92);  % in original 16-bit value
%threshold = double(600)/65535;  % between 0-1
threshold = 0.02;

% binary image
Ibw2 = im2bw(I, threshold); % MAKING THE IMAGE BINARY
Ibw2 = bwareafilt(Ibw2,[1 25]); % function to set the pixel size of the RCPs to include in the 
figure, imshow(Ibw2, []);

% find peaks
figure()
Idil = imdilate(I, ones(5)); % IMAGE DIALATION
centroid = I == Idil & Ibw2;
[Y,X] = ind2sub(size(I), find(centroid));
figure; subplot(121);
imshow(I,[]); 
hold on; plot(X, Y, 'r.');

% keep only isolated
idx = rangesearch([X,Y], [X,Y], maxdist_between_spots);
nNN = cellfun(@length, idx); % apply function to each cell in cell array
plot(X(nNN==1), Y(nNN==1), 'yo');

subplot(122); hold on;
bell_dist = 10;
fid = fopen(filename, 'w'); % fopen is used to open file or obtain information about open files
% the w will specify the write access to the input file 
for i = find(nNN==1)'
    try
        plot(1:(2*bell_dist+1), I(Y(i),X(i)-bell_dist:X(i)+bell_dist) - min(I(Y(i),X(i)-bell_dist:X(i)+bell_dist)), 'k-');
        drawnow;
        fprintf(fid, lineformat('%d', 2*bell_dist+1), I(Y(i),X(i)-bell_dist:X(i)+bell_dist));
    end
end
fclose(fid);


function fmt = lineformat(vartype, repeat)
% fmt = lineformat(vartype, repeat)
% line format for csv file writing
% ends with newline character
% Xiaoyan, 2017


fmt = repmat([vartype, ','], 1, repeat);
fmt = [fmt(1:end-1), '\n'];

end
