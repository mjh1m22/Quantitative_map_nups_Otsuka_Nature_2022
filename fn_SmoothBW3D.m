function [outStack] = fn_SmoothBW3D(inStack, len)
inTotPixel = sum(sum(sum(inStack)));
inStack = smooth3(inStack, 'box', len);
outStack = double(inStack>=0.5);
th = 0.49;
while sum(sum(sum(outStack)))<inTotPixel*0.99 && th>0.2
    outStack = double(inStack>th);
    th = th-0.01;
end
end

