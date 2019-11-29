function rgb = rgbscale(xyz)
% rgb = rgbscale(xyz)
% Xiaoyan, 2018

rgb = xyz - min(xyz, [], 1);
rgb = bsxfun(@rdivide, rgb, max(rgb, [], 1));
end
