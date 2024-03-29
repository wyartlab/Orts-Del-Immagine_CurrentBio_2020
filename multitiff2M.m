function M = multitiff2M(tiff,framesel);
tiffinfo = imfinfo(tiff);
%framesel= total number of acquired frames;

numframes = length(tiffinfo);
%nb de frames

M = zeros(tiffinfo(1).Height,tiffinfo(1).Width,['uint' num2str(tiffinfo(1).BitDepth)]);

for frame=framesel;
    curframe = im2uint8(imread(tiff,frame));
    M(:,:,frame-framesel(1)+1) = curframe;
   
end
