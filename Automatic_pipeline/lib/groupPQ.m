function Group = groupPQ(Idx_PQ)
% for droplet barcoding
% recursive grouping
% Xiaoyan, 2016



% GROUP reads (clustering)
Group = zeros(size(Idx_PQ,1),1,'double');

c = 0;
while ~isempty(Idx_PQ)
    c = c+1;
    FindP(Idx_PQ(1,2));
end


% child functions
    function FindP(select_P)
        bin_P = ismember(Idx_PQ(:,2),select_P);
        
        if nnz(bin_P)
            
            select_Q = Idx_PQ(bin_P,3);
            
            Group(Idx_PQ(bin_P,1)) = c;
            Idx_PQ(bin_P,:) = [];
            
            FindQ(unique(select_Q));
        end
        
    end

    function FindQ(select_Q)
        bin_Q = ismember(Idx_PQ(:,3),select_Q);
        if nnz(bin_Q)
            select_P = Idx_PQ(bin_Q,2);
            
            Group(Idx_PQ(bin_Q,1)) = c;
            Idx_PQ(bin_Q,:) = [];
            
            FindP(unique(select_P));
        end
        
    end
end

