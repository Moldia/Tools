function new_pos = correctcoordinates_f(pos,resize_factor)
% correct coordinates for plotting and image processing


new_pos = (pos-1/resize_factor/2+.5)*resize_factor+1;

end
