green = "#BCBF65";
blue = "#51A6A6";
muted_red = "#D37D7D"; % Muted red color

figure(1)
set(gcf, 'Color', 'w'); % Set background color to white

tiledlayout(2,3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Adjust spacing

% Store axis handles
ax = gobjects(3,3); % 3 rows, 3 columns

% Function to plot with horizontal reference lines
function ax_handle = plot_with_lines(x, y, struct, xlabel_text, color)
    muted_red = "#D37D7D"; % Muted red color
    ax_handle = gca; % Get current axes
    plot(x, y, 'color', color, 'LineWidth', 1.5);
    hold on;
    yline(struct.L, 'color', muted_red, 'LineStyle', ':', 'LineWidth', 1.2);
    yline(struct.L-struct.Hlake, 'color', muted_red, 'LineStyle', ':', 'LineWidth', 1.2);
    yline(0, 'color', muted_red, 'LineStyle', ':', 'LineWidth', 1.2);
    
    % Add text labels next to horizontal lines
    text(min(x), struct.L, sprintf('n= %.3f', struct.n1), 'Color', muted_red, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(min(x), struct.L-struct.Hlake, sprintf('n= %.3f', struct.n2), 'Color', muted_red, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(min(x), 0, sprintf('n= %.3f', struct.n3), 'Color', muted_red, 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    
    hold off;
    grid on;
    xlabel(xlabel_text, "FontSize", 12);
end

% First row: n_{gas} vs. height
nexttile;
ax(1,1) = plot_with_lines(a.n, a.z, a, 'n_{gas}', blue);
ylabel('height in column m', "FontSize", 12);

nexttile;
ax(1,2) = plot_with_lines(b.n, b.z, b, 'n_{gas}', blue);

nexttile;
ax(1,3) = plot_with_lines(c.n, c.z, c, 'n_{gas}', blue);

% Second row: rho vs. height
nexttile;
ax(2,1) = plot_with_lines(a.rho, a.z, a, 'mixture density kg/m3', green);
ylabel('height in column m', "FontSize", 12);

nexttile;
ax(2,2) = plot_with_lines(b.rho, b.z, b, 'c m/s', green);

nexttile;
ax(2,3) = plot_with_lines(c.rho, c.z, c, 'c m/s', green);






% Third row: c vs. height
nexttile;
ax(3,1) = plot_with_lines(a.c, a.z, a, 'c m/s', green);
ylabel('height in column m', "FontSize", 12);

nexttile;
ax(3,2) = plot_with_lines(b.c, b.z, b, 'c m/s', green);

nexttile;
ax(3,3) = plot_with_lines(c.c, c.z, c, 'c m/s', green);

% Link y-axes so they stay aligned
linkaxes(ax, 'y');

% Remove y-axis tick labels from non-leftmost plots
set(ax(:,2:3), 'YTickLabel', []);
