function [] = vis_pfm_diff(pfmAPath, pfmBPath)
%VIS_PFM_DIFF Visualize percent difference between PDMs

	FONT_SIZE = 24;

	a = parsePfm(pfmAPath);
    b = parsePfm(pfmBPath);

    diff = abs(a - b);
    maxAbsVal = max(  max(abs(a(:))), max(abs(b(:))) );
    
    imshow((diff / maxAbsVal) * 100.0);
    colorbar;
    colormap jet;
    t = title('Absolute percent diff');
    

end

