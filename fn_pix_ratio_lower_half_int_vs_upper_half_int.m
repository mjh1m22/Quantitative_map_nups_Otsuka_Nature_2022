function [pix_ratio_lhi_vs_uhi] = fn_pix_ratio_lower_half_int_vs_upper_half_int(nucRim, nup_stack, ratio_thresh)
%% This function returns the ratio of pixels containing lower half by the upper half of total intensity on nuclear surface 

[xSize, ySize, zSize] = size(nup_stack);

maxVal = max(max(max(nup_stack)));
histN = round(maxVal)+1;
histogram = zeros(1,histN);

for i = 1:xSize
    for j = 1:ySize
        for k = 1:zSize
            if nucRim(i,j,k) >0
                hIndex = round(nup_stack(i,j,k))+1;
                histogram(hIndex) = histogram(hIndex)+1;
            end
        end
    end
end

tot_int = sum(histogram .*[1:histN]);
tot_pix = sum(histogram);

for i = histN:-1:1
    upper_tot_int = sum(histogram(i:histN) .*[i:histN]);
    if upper_tot_int/(tot_int-upper_tot_int) >= ratio_thresh
        upper_tot_pix = sum(histogram(i:histN));

        break;
    end    
end
pix_ratio_lhi_vs_uhi = (tot_pix - upper_tot_pix)/upper_tot_pix;
end