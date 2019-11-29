function [name, varargout] = removereads(name, nameExclude, varargin)
% [name, varargout] = removereads(name, nameExclude, varargin)
% remove reads
% Xiaoyan, 2017


if ~iscell(nameExclude)
    nameExclude = {nameExclude};
end
[uniName, ~, idxName] = unique(name);

idxNameExclude = [];
for i = 1:length(nameExclude)
    idxNameExclude = [idxNameExclude;...
        find(cellfun(@(v) strcmp(v, nameExclude{i}), uniName))];
end

readsExclude = ismember(idxName, idxNameExclude);

name = name(~readsExclude);
for i = 1:length(varargin)
    varargout{i} = varargin{i}(~readsExclude,:);
end

end
