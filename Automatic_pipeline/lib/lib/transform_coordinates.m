function new_pos = transform_coordinates(pos, scaling, imsize, ang, yup, xleft,...
    imsize_final, imscale)
% new_pos = transform_coordinates(pos, scaling, imsize, ang, yup, xleft,...
%     imsize_final, imscale)
% Xiaoyan, 2018

% scaling first
pos_re = pos*scaling;

% center as origo
pos_transformed = bsxfun(@minus, pos_re(:,1:2), fliplr(imsize/imscale)/2);

% rotation
rotation_angle = -1*ang/180*pi;
rototation_mat = [...
    cos(rotation_angle), sin(rotation_angle);...
    -sin(rotation_angle), cos(rotation_angle)];
new_pos = pos_transformed*rototation_mat;

% origo
new_pos = bsxfun(@plus, fliplr(imsize_final/imscale)/2, new_pos);

% translation
new_pos(:,1) = new_pos(:,1)-xleft/imscale;
new_pos(:,2) = new_pos(:,2)-yup/imscale;

end

