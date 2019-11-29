% rename and move reference images to AlignedImages folder


ref = 'base1';


mkdir('AlignedImages');

for c = 1:6
    movefile([ref, '\img_t', num2str(c), '_z1_c1.tif'],...
        ['AlignedImages\', ref, '_c', num2str(c), '_ORG.tif']);
end
