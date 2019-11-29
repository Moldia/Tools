function newobjprop = renumberedgeobj...
    (ntiles, correctiontable, objlist, column)
% newobjprop = renumberedgeobj...
%     (ntiles, correctiontable, objlist, column)
% renumber ojects at the edge based on the correction table
% Xiaoyan, 2017


newobjprop = zeros(size(objlist,1), 1);
for i = 1:ntiles(2)
    for j = 1:ntiles(1)        
        t = ntiles(1)*(i-1)+j;
        from = objlist(objlist(:,3)==t, column);
        if ~isempty(from)
            [uniFrom, ~, idxFrom] = unique(from);
            to = correctiontable{i,j};
            for k = 1:length(uniFrom)
                if uniFrom(k)
                    try
                        from(idxFrom==k) = to(uniFrom(k), 2);
                    catch
                        from(idxFrom==k) = -1;   % cell not found in the object image
                    end
                end
            end
            newobjprop(objlist(:,3)==t) = from;
        end
    end
end

end
