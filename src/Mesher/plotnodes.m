%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose
% =======
% Plot nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function plotnodes(xe, coord) 
 

close all;
x = coord(:,1);

width = 6;     % Width in inches
height = 4;    % Height in inches
alw = 5.0;     % AxesLineWidth
fsz = 18;      % Fontsize
lw = 2.0;      % LineWidth

figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
axis([min(xe) max(xe) -1.0/5.0 1.0/5.0 ]);
xlabel('x','FontSize', fsz);
        
% 
plot(xe, zeros(length(xe),1), 'mo', 'LineWidth', lw, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.49 1 .63], ...
    'MarkerSize', 15);

plot(x, zeros(length(x),1),'x','Linewidth', lw, ...
    'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',15);

legend('Element Boundary Nodes', 'FE Nodes')

return
end
