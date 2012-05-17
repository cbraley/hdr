function [] = plot_ctf(ctfFile, pointsFile)
%PLOT_CTF Plot CTF and (optionally) the points used to generate it.
%   ctfFile is the path to the CTF data.
%   pointsFile is the path to the points data.  Use '' to ignore this
%   parameter.

    %Constants
    LINE_W    = 3;
    LINE_COLOR = [0 0 0];
    DOT_SIZE  = 5;
    DOT_COLOR = [0 0 1];
    FONT_SIZE = 24;
    PERC = .1;

	%Maybe plot the points
    hold on;
    if strcmp(pointsFile, '') ~= 1
        points = load(pointsFile);
        plot(points(:,1), points(:,2), 'o', 'MarkerEdgeColor', DOT_COLOR, ...
            'MarkerSize', DOT_SIZE);
    end
    
    %Plot the CTF
    pix_vals = 0:255;
    ctf = load(ctfFile);
    plot(pix_vals, ctf, '-', 'LineWidth', LINE_W, 'Color', LINE_COLOR);
	

    %Make the plot pretty
    t = title(sprintf('CTF From file: %s', ctfFile));
    set(t, 'FontSize', FONT_SIZE);
    axis tight;
    grid on;
    grid minor;
    t = xlabel('8 bit pixel values');
    set(t, 'FontSize', FONT_SIZE);
    t = ylabel('Exposure');
    set(t, 'FontSize', FONT_SIZE);
    hold off;
    
    %Make the plot fit the CTF well
    ctfMin = min(ctf);
    ctfMax = max(ctf);
    ctfRange = ctfMax - ctfMin;
    ylim([ctfMin - (PERC * ctfRange), ctfMax + (PERC * ctfRange)]);
    
    
end

