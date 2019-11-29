function K = ripleysk_modified(locs,h,box,method,ROIbw)
% RipleysK: Calculate K statistic
% 
% K = RipleysK(locs,dist, box,method) calculates G, the K statistic at each 
% distance for the data with x-y coordinates given by locs, and the
% bounding rectangle given by box=[minX maxX minY maxY].
% If method=0, no edge correction is applied.
% If method=1, points are only used if they are at least h units away from
% the edge.
% If method=2, an ROI defineded by ROI is used instead a square, with eget
% correction.
%
% Note: The L statistic may be calculated from the K statistic as follows: 
%   L = sqrt(K/pi)-h;
%   
% modified by Xiaoyan

if nargin<4, method=1; end
[N,k] = size(locs);
if k~=2, error('locs must have two columns'); end
rbox = min([locs(:,1)'-box(1);box(2)-locs(:,1)';locs(:,2)'-box(3); box(4)-locs(:,2)'] );
% rbox is distance to box

DX = repmat(locs(:,1),1,N)-repmat(locs(:,1)',N,1);
DY = repmat(locs(:,2),1,N)-repmat(locs(:,2)',N,1);
DIST = sqrt(DX.^2+DY.^2);
DIST = sort(DIST);

if method==1
    K = zeros(length(h),1);
    for k=1:length(K)
        I = find(rbox>h(k));
        
        if ~isempty(I)
            K(k) = sum(sum(DIST(2:end,I)<h(k)))/length(I);
        else
            K(k) = -1;
        end
        
    end
elseif method==2
    K = zeros(length(h),1);    
    ROIdist = bwdist(~ROIbw);
    for k=1:length(K)
        Irec = find(rbox>h(k));
        I = [];
        for i = 1:length(locs)
            if ROIdist(floor(locs(i,2))+1,floor(locs(i,1))+1)>h(k)
                I = [I;i];
            end
        end        
        I = I(ismember(I,Irec));
        
        if ~isempty(I)
            K(k) = sum(sum(DIST(2:end,I)<h(k)))/length(I);
        else
            K(k) = -1;
        end
    end
elseif method==0
    K = zeros(length(h),1);
    for k=1:length(K)
        K(k) = sum(sum(DIST(2:end,:)<h(k)))/N;
    end
end

if method == 2
    ROI_range = ROIbw(floor(box(3))+1:floor(box(4))+1,floor(box(1))+1:floor(box(2))+1);
    lambda = N/nnz(ROI_range);
else
    lambda = N/((box(2)-box(1))*(box(4)-box(3)));
end
   
K = K/lambda;