function [img, scaleFactor] = parsePfm(filePath)
%parsePfm Parses .pfm images.
%   Bring a .pfm image into matlab.

    %Open the file and check for errors
    fid = fopen(filePath, 'r');
    if fid == -1
        error('pfm:IOError', 'Could not open file!');
    end
    
    %File opened OK
    %
    %.pfm headers have 3 ASCII lines
    %Line 1: The text 'PF' or 'Pf' where the former denotes color and the
    % latter denotes grayscale
    %Line 2: Two integers, width then height.
    %Line 3: A single signed decimal number S
    %    is S < 0 then the file is little endian
    %    otherwise the file is big endian
    %    |S| is a scale factor to relate pixel samples to a physical
    %    quantity(like radiance for example).  
    %
    
    %Info to determine during header parse
    numChannels = 0; %1 = grayscale, 3 = RGB
    imWidth     = 0;
    imHeight    = 0; 
    isBigEndian = 0;
    scaleFactor = 1; %Described above, ignored by this code for now
    
    %Read the whole 3 line header
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    line3 = fgetl(fid);
    if ~(ischar(line1) && ischar(line2) && ischar(line3))
        fclose(fid);
        error('pfm:IOError', 'Header was incomplete!');
    end
    
    %Parse line 1, determine color or BW 
    if strcmp(line1, 'PF') == 1 %Color
        numChannels = 3;
    elseif strcmp(line1, 'Pf') == 1 %Gray
        numChannels = 1;
    else %Invalid header
        fclose(fid);
        error('pfm:IOError', 'Invalid .pfm header!');
    end
    
    %Parse line 2, get image dims
    [dims, foundCount, errMsg] = sscanf(line2, '%u %u');
    if numel(dims) ~= 2 || strcmp(errMsg,'') ~= 1 || foundCount ~= 2
        fclose(fid);
        error('pfm:IOError', 'Dimensions line was malformed!');
    end
    imWidth  = dims(1);
    imHeight = dims(2);
    
    %Line 3, the endianness+scale line
    [scale, matchCount, errMsg] = sscanf(line3, '%f');
    if matchCount ~= 1 || strcmp(errMsg,'') ~= 1
        fclose(fid);
        error('pfm:IOError', 'Endianness+Scale line was malformed!');
    end
    scaleFactor = abs(scale);
    endianChar = 'n';
    if scale < 0.0
        isBigEndian = 0;
        endianChar = 'l';
    else
        isBigEndian = 1;
        endianChar = 'b';
    end
    
    %Allocate image buffer
    img = zeros(imWidth, imHeight, numChannels);
    totElems = numel(img);
    
    %Now at last parse in the pixel raster
    %the raster is a 4 byte valeues arranged left to right, starting in the
    %upper left corner of the image.  In the case of a color image,
    %channels are interleaved 
    [rawData, numFloatsRead] = fread(fid, totElems, 'single', 0, endianChar);
    if numFloatsRead ~= totElems
        fclose(fid);
        error('pfm:IOError', 'Raster data did not match header!');
    end
    fclose(fid);
    
    %Put the data into the output buffer after re-shaping to match Matlab's
    %ordering and color channel order.
    imDataInReadOrder = reshape(rawData, [numChannels, imWidth, imHeight]);
    img = permute(imDataInReadOrder, [3, 2, 1]);
    
end
