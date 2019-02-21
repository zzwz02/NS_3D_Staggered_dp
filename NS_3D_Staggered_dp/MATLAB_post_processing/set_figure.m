function [gca, gcf] = set_figure(gca, gcf, width, height, font_size)

set(gca,'FontSize',font_size)
set(gcf, 'Units', 'centimeters', 'PaperUnits', 'centimeters', 'Color', [1 1 1]);
temp=get(gcf, 'Position');
set(gcf, 'Position', [temp(1), temp(2), width, height], 'PaperSize', [width, height])

end

