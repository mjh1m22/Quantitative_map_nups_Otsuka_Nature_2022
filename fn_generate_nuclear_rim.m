function [nucRim] = fn_generate_nuclear_rim(nucRegion, numDil, numErd, strElRim)
nucRimDil = nucRegion;
for dilIdx = 1:numDil
    nucRimDil = imdilate(nucRimDil,strElRim);
end
nucRimErd = nucRegion;
for erdIdx = 1:numErd
    nucRimErd = imerode(nucRimErd,strElRim);
end
nucRim = double(nucRimDil -nucRimErd);
end