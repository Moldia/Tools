function [flo_rotate, Ifuse_rotate] = rotateimage(flo,ang,ref)

flo_rotate = imrotate(flo,ang);
Ifuse_rotate = imfuse(flo_rotate,ref);

end
