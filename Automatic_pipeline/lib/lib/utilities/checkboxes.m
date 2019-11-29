function cbvalues = checkboxes(names)
% cbvalues = checkboxes(names)
% create checkboxes and get user selection
% Xiaoyan, 2017

names = unique(names, 'stable');
cbh = zeros(99, ceil(numel(names)/99)); % 99 checkboxes per figure window
pbh = zeros(1, ceil(numel(names)/99));  % one pushbutton per figure window
cbvalues = []; 

fnum = 1;
f(fnum) = figure('visible', 'off');

for k = 1:numel(names)
    if ceil(k/99) > fnum
        
        fnum = fnum + 1;
        f(fnum) = figure('visible', 'off');
    end
    
    x = ceil(k/20);
    y = k - 20*(x-1);
    x = x - 5*(fnum-1);
    
    % checkboxes
    cbh(mod(k,99)+1, fnum) = uicontrol(f(fnum), 'Style', 'Checkbox',...
        'String', names{k},...
        'Value', 1,...
        'Units', 'Normalized', 'Position', [.2*(x-1) 1-.05*y .2 .05],...
        'Callback', {@cb_callback, k});
end

% make figures visible
for i = fnum:-1:1
    set(f(i), 'visible', 'on');
end

for i = 1:fnum
    figure(f(i));
    % pushbutton
    pbh(i) = uicontrol(f(i), 'Style', 'pushbutton', 'String', 'Confirm',...
        'Units', 'Normalized', 'Position', [.82 0 .16 .05],...
        'Callback', @get_cbvalues);    
    uiwait(f(i));
end

    function cb_callback(hObject, eventData, checkBoxId)
        value = get(hObject, 'Value');
    end


    function get_cbvalues(hObject, eventData)
        % get all checkbox values from current figure
        ch = get(gcf, 'child');
        for n = numel(ch):-1:1
            if strcmp(ch(n).Style, 'checkbox')
                cbvalues = [cbvalues; ch(n).String, num2cell(ch(n).Value)];
            end
        end
        close(gcf);
        drawnow;
    end

end

