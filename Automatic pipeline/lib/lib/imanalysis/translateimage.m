function Ifuse_translate = translateimage(yup, xleft, flo_rotate, ref)
% Xiaoyan, 2017

if yup>=0
    if xleft>=0
        Ifuse_translate = imfuse(flo_rotate(yup+1:end,xleft+1:end), ref);
    else
        Ifuse_translate = imfuse(flo_rotate(yup+1:end,:), ref(:,-xleft+1:end));
    end
else
    if xleft>=0
        Ifuse_translate = imfuse(flo_rotate(:,xleft+1:end), ref(-yup+1:end,:));
    else
        Ifuse_translate = imfuse(flo_rotate, ref(-yup+1:end,-xleft+1:end));
    end
end
