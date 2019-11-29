function seqdata = decodinginput(inputfile, colChannels, colOthers)
disp('Start decoding input')
% seqdata = decodinginput(inputfile, colChannels, colOthers)
% return sequencing data matrix from CP output
% Xiaoyan 2014-11-28

% input check 
temp = 3 + length(colChannels) + length(colOthers);
if temp ~= 12
    error('Wrong column number detected.');
end

% import
seq = importdata(inputfile, ',', 1);


blobs = size(seq.data, 1);

% preallocate
seqdata = zeros(blobs, 14);

% check column numbers
column1_exist = find(colChannels);
column1_missing =  find(colChannels==0);
column2_exist = find(colOthers);
column2_missing = find(colOthers==0);

if colOthers(1)==0 || colOthers(2)==0
    error('X and Y position mush exist.')
end
% reform into output matrix
opts = detectImportOptions(inputfile);
table=readtable(inputfile,opts);

seqd=(seq.data);
a=[seqd(:,1),seqd(:,2)];
a=array2table(a);
seq.data=horzcat(table(:,1:11),a);
seq.data=table2array(seq.data);
table2=table2array(table);
seqdata(:,1:11) = seq.data(:,1:11);
%seqdata(:,12:13) = table2(:,7:8);
seqdata(:,12) = table2(:,8);
seqdata(:,13) = table2(:,7);
%

seqdata(:,9) = zeros(blobs,1);
ssd=size(seqdata);
seqdata(:,14) = [1:ssd(1)];

seqdata(:,column1_exist-1) = ...
    seq.data(:,colChannels(column1_exist));
disp('Finish decoding input')
end
            
        
  
  

