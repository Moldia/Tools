function blobinpolygon = binreads(Pos, grid_size)
% blobinpolygon = binreads(Pos, grid_size)
% group reads into square
% Xiaoyan, 2018

% count reads in grid
grid_nrx =  ceil(max(Pos(:,1))/grid_size);
grid_nry =  ceil(max(Pos(:,2))/grid_size);
grid_nr = grid_nrx*grid_nry;

% scaling
Pos_pixel = round(correctcoord(Pos,1/grid_size));

blobinpolygon = false(size(Pos,1),grid_nr);
Counted = false(length(Pos),1);

for j = 1:grid_nrx   % j along x axis, column
%     fprintf('%.2f%s\n',double((j-1)/grid_nrx)*100,'% grid processed.');
    
    for i = 1:grid_nry  % i along y axis, row
        temp_in = Pos_pixel(:,1)==j & Pos_pixel(:,2)==i;
        temp_in = temp_in & ~Counted;
        
        % blobs within grid (logical)
        k = (j-1)*grid_nry+i; % counting direction: y
        blobinpolygon(:,k) = temp_in;
        
        % already counted blobs
        Counted = Counted | temp_in;
        
    end
end

disp('100.00% grid processed.');
end
