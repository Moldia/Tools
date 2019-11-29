function [CellCorrectionTable, CellProps, imOutlines] = correct_edge_objects...
    (imLabel, imOutlines, cellprop, tilepos, Ivis, resizef)
% [CellCorrectionTable, CellProps, imOutlines] = correct_edge_objects...
%     (imLabel, imOutlines, cellprop, tilepos, Ivis, resizef)
% correct edge objects
% Xiaoyan, 2017

tilesize = size(imLabel{1},1);
ntilesY = size(imLabel,1);
ntilesX = size(imLabel,2);
startlabel = 0;

iempty = zeros(tilesize, tilesize, 'uint16');
CellCorrectionTable = cell(ntilesY,ntilesX);
CellProps = [];

% visualization
figure, imshow(Ivis)
hold on

for i = 1:ntilesY
    fprintf('%.2f percet finished.\n', (i-1)/ntilesY*100);
    for j = 1:ntilesX
        t = ntilesX*(i-1)+j;
       
        I = imLabel{i,j};
                    
        try
            I1 = imLabel{i,j-1};
        catch
            I1 = iempty;
        end
        
        try
            I2 = imLabel{i-1,j};
        catch
            I2 = iempty;
        end
    
        % find cut objects - x
        temp = I1(:,end)~=0 & I(:,1)~=0;
        % outline images
        if nnz(temp)
            imOutlines{i,j-1}(temp,end,:) = 0;
            imOutlines{i,j}(temp,1,:) = 0;
        end
        % cut objects
        temp = [I1(temp,end),I(temp,1)];
        matchx = unique(temp, 'rows');
        
        % find cut objects - y
        temp = I2(end,:)~=0 & I(1,:)~=0;
        % outline images
        if nnz(temp)
            imOutlines{i-1,j}(end,temp,:) = 0;
            imOutlines{i,j}(1,temp,:) = 0;
        end

        temp = [I2(end,temp); I(1,temp)]';
        matchy = unique(temp, 'rows'); 

        % corrected cell properties
        tileCellProps = correctprops(t, I, cellprop, tilepos);
        CellProps = [CellProps; tileCellProps];

        % renumber objects
        I = uint32(I);  % in case total object number exceed 16-bit      
        I(I~=0) = I(I~=0) + startlabel;
        uniI = unique(I(:));
        CellCorrectionTable{i,j} = [(1:length(uniI(2:end)))', uniI(2:end)];
        
        if max(I(:))
            startlabel = max(I(:));
        end

        % correct cut objects
        if ~isempty(matchx) 
            % visualization
            temp = CellProps(CellProps(:,3)==t & ismember(CellProps(:,2),matchx(:,2)),end-1:end);
            plot(temp(:,1)*resizef,temp(:,2)*resizef,'<');
            temp = CellProps(CellProps(:,3)==t-1 & ismember(CellProps(:,2),matchx(:,1)),end-1:end);
            plot(temp(:,1)*resizef,temp(:,2)*resizef,'>');
            
            from = CellCorrectionTable{i,j-1};
            to = CellCorrectionTable{i,j};
            
            unifrom = unique(matchx(:,1));
            for k = 1:length(unifrom)
                temp = matchx(matchx(:,1)==unifrom(k),2);
                to(temp,2) = to(temp(end),2);
            end            
            
            for k = 1:size(matchx,1)
                from(from(:,1)==matchx(k,1),2) = to(matchx(k,2),2);
            end
            
            CellCorrectionTable{i,j} = to;
            CellCorrectionTable{i,j-1} = from;
        end
        
        if ~isempty(matchy)
            % visualization
            temp = CellProps(CellProps(:,3)==t & ismember(CellProps(:,2),matchy(:,2)),end-1:end);
            plot(temp(:,1)*resizef,temp(:,2)*resizef,'^');           
            temp = CellProps(CellProps(:,3)==t-ntilesX & ismember(CellProps(:,2),matchy(:,1)),end-1:end);
            plot(temp(:,1)*resizef,temp(:,2)*resizef,'v');
            
            from = CellCorrectionTable{i-1,j};
            to = CellCorrectionTable{i,j};
            
            unifrom = unique(matchy(:,1));
            for k = 1:length(unifrom)
                temp = matchy(matchy(:,1)==unifrom(k),2);
                to(temp,2) = to(temp(end),2);
            end
            
            for k = 1:size(matchy,1)
                from(from(:,1)==matchy(k,1),2) = to(matchy(k,2),2);
            end
            
            CellCorrectionTable{i,j} = to;
            CellCorrectionTable{i-1,j} = from;
        end
    end
end
fprintf('100 percet finished.\n');

clear I temp iempty
end

function reCellProps = correctprops(t, I, cellprop, tilepos)
uniI = unique(I(:));
if nnz(cellprop(:,3)==t) ~= length(uniI)-1
    warning(['tile', num2str(t), ' object numbers do not match: ',...
        num2str(nnz(cellprop(:,3)==t)), ' vs ', num2str(length(uniI)-1)])
    temp = regionprops(I, 'Area', 'Centroid');
    centroid = cat(1,temp.Centroid)+repmat(tilepos(t,2:3), size(temp,1),1);
    % THIS PART NOT FINISHED, NOT TESTED
    [nn, D] = knnsearch(centroid, cellprop(cellprop(:,3)==t,end-1:end), 'K', 1);
    temp = [cat(1,temp.Area).*double(D<=50), centroid.*repmat(double(D<=50),1,2)];
    reCellProps = [repmat(t,size(temp,1),1), nn.*double(D<=50),...
        repmat(t,size(temp,1),1),temp];
else
    reCellProps = cellprop(cellprop(:,3)==t,:);
end

end

