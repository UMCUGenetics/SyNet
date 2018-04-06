function [box_h, patch_h] = BoxPlotEx(Dist, varargin)
box_h = boxplot(Dist, varargin{:});
set(box_h, 'LineWidth', 2);
XData = get(box_h, 'XData');
YData = get(box_h, 'YData');
Box_Color = get(box_h, 'Color');
patch_h = patch(XData{6}([1 2 2 1]), YData{5}([1 1 2 2]), 'b', 'FaceColor', Box_Color{5}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
median_lines = findobj(box_h, 'type', 'line', 'Tag', 'Median');
set(median_lines, 'Color', Box_Color{5}*0.8);
set(findobj(box_h,'tag','Outliers'), 'MarkerEdgeColor', Box_Color{5}*0.7)
% uistack(box_h,'top');
end