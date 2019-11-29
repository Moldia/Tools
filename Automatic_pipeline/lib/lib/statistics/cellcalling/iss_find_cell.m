function cells = iss_find_cell(cellpos, o)
% cells = iss_find_cell(cellpos, o)

cells = knnsearch(o.CellYX, fliplr(cellpos), 'K', 3);

end