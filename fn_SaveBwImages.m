function fn_SaveBwImages(imStack, outDir, fnPre)
%Saves bwImages for amira rendering
Nz = size(imStack,3);
if ~exist(outDir)
    mkdir(outDir);
end

dumImg = zeros(size(imStack(:,:,1)));
for f = 0:Nz+1
    saveFileName = [outDir fnPre num2str(f, '%03.f') '.tif'];
    if f == 0 || f == Nz+1
        imwrite(dumImg, saveFileName, 'tif');
    else
        imwrite((imStack(:,:,f)), saveFileName, 'tif');
    end
end
end