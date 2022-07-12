function [thresh3D, thresh2D] = fn_PC_GetThresholds(stack)
% This function calculate 3D and 2D threshold of the stack
Nz = size(stack, 3);
thresh2D = zeros(Nz,1);
divFact = max(max(max(stack)))/1000;
if divFact > 1
    stackDiv = stack/divFact;
    [thresh3D, histChr] = Otsu_3D_Img(stackDiv, 0);
    thresh3D = thresh3D * divFact;
else
    stackDiv = stack;
    [thresh3D, histChr] = Otsu_3D_Img(stackDiv, 0);
end

%Update threshold for different slices
for i = 1:Nz
    thresh2D(i) = Otsu_2D(stackDiv(:,:,i));
    if divFact > 1
        thresh2D(i) = thresh2D(i) * divFact;
    end
end
thresh2D = smooth(thresh2D, 15);
end

