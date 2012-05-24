function [] = vis_HDR(hdrPfmPath, Npath, residualPath )
%VISUALZE_HDR_RESULTS Visaulize output of gen_hdr.
%   hdrPfmPath is a path to a HDR .pfm image. May be '' to skip.
%   Npath is a path to a LDR image generated via --out_n . May be '' to skip.
%   residualPath is a path to a LDR image generated via --out_r . May be '' to skip.

    FONT_SIZE = 24;

    %Potentially show the monochrome HDR
    if( strcmp(hdrPfmPath,'') ~= 1)
        hdr = parsePfm(hdrPfmPath);
        figure;
        
        %Check for RGB (unsupported)
        if( size(hdr,3) ~= 1)
            fprintf('Error - Only monochrome HDRs are currently supported!');
            return;
        end
        
        imshow(hdr);
        colormap jet;
        colorbar;
        axis image;
        t = title(sprintf('HDR image from %s',hdrPfmPath));
        set(t, 'FontSize', FONT_SIZE);
    end
    
    %Potentialyl show # samples
    if( strcmp(Npath,'') ~= 1)
        im = imread(Npath);
        figure;
        image(im,'CDataMapping','scaled');
        temp = colormap('jet');
        colormap(flipud(temp));
        colorbar;
        axis image;
        t = title(sprintf('Num Valid HDR samples visualization from %s',Npath));
        set(t, 'FontSize', FONT_SIZE);
    end
    
	%Plot residual
    if( strcmp(residualPath,'') ~= 1)
        im = parsePfm(residualPath);
        figure;
        image(im,'CDataMapping','scaled');
        colormap jet;
        colorbar;
        axis image;
        t = title(sprintf('Residual visualization from: %s',residualPath));
        set(t, 'FontSize', FONT_SIZE);
    end


end

