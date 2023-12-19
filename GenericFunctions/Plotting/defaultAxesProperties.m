function defaultAxesProperties(ax, offsetFlag)

set(ax, 'TickDir','out', 'TickLength', [0.01 0.001]...
    ,  'color', 'none', 'box','off', 'XColor', 'k', 'YColor', 'k',...
    'FontName', 'Calibri', 'LineWidth', 0.5)

if nargin>1
    if offsetFlag
        offsetAxes(ax);
    end
end