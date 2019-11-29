function insituQC(name, beforeQT, afterQT, nHybs, fontsize)
% insituQC(name, beforeQT, afterQT, nHybs, fontsize)
% in situ seq quality check
% output image
% Xiaoyan, 2016-2-18

warning off MATLAB:MKDIR:DirectoryExists
mkdir('QC');

num_bases = 4;
[Y{nHybs:-1:1}] = ndgrid(1:num_bases);
nummatrix = reshape(cat(nHybs+1, Y{:}), [], nHybs);
numlist = zeros(size(nummatrix,1), 1);
for i = 1:nHybs
    numlist = numlist + nummatrix(:,i)*10^(nHybs-i);
end
list = num2barcode(numlist);
clear Y

Reads = zeros(length(list),2);

% before
count = importdata(beforeQT, ',');
count(1,:) = [];
count = cellfun(@(v) strsplit(v, ','), count, 'uni', 0);
reads = cellfun(@(v) v{:,1}, count, 'uni', 0);
counts = cellfun(@(v) v{:,2}, count, 'uni', 0);
counts = cell2mat(cellfun(@str2double, counts, 'uni', 0));
names = cellfun(@(v) v{:,3}, count, 'uni', 0);

idx_NNNN = strcmp(names, 'NNNN');
Expected_reads = reads(~idx_NNNN);
Expected_names = names(~idx_NNNN);

Expected = zeros(nnz(~idx_NNNN), 2);
Expected(:,1) = counts(~idx_NNNN);

existreads = ismember(list,reads);
Reads(existreads,1) = counts;

% after
count = importdata(afterQT,',');
count(1,:) = [];
count = cellfun(@(v) strsplit(v,','),count, 'uni', 0);
reads = cellfun(@(v) v{:,1},count, 'uni', 0);
counts = cellfun(@(v) v{:,2},count, 'uni', 0);
counts = cell2mat(cellfun(@str2double,counts, 'uni', 0));
names = cellfun(@(v) v{:,3},count, 'uni', 0);

idx_NNNN = strcmp(names,'NNNN');
idx_expected = ismember(Expected_reads,reads(~idx_NNNN));
Expected(idx_expected,2) = counts(~idx_NNNN);

existreads = ismember(list,reads);
Reads(existreads,2) = counts;

clf; 
subplot(1,2,1),hold on;
plot(Reads(:,1),Reads(:,2),'.')
plot(Expected(:,1),Expected(:,2),'o');
text(Reads(:,1),Reads(:,2),list,'fontsize',fontsize);
set(gca,'YScale','log','XScale','log');
xlabel('before QT');
ylabel('after QT');
pbaspect([1,1,1]);

name = strrep(name,'_','\_');
title(name);

subplot(1,2,2),hold on;
plot(Reads(:,1),Reads(:,2),'.')
plot(Expected(:,1),Expected(:,2),'o');
text(Expected(:,1),Expected(:,2),Expected_names,'fontsize',fontsize);
set(gca,'YScale','log','XScale','log');
xlabel('before QT');
ylabel('after QT');
pbaspect([1,1,1]);

title(name);

name = strrep(name,'\_','_');
print(gcf,'-dpng','-r1000',['QC\' name '.png']);

end
