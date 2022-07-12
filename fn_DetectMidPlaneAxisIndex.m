function [refAxIdNuc, axId1, axId2] = fn_DetectMidPlaneAxisIndex(refAxIdNuc, V)
minAng = 90;
minAngIndex = 1; %Otherwise chose the one that has minimum angle with previous
for comp = 1:3
    ang = acos(abs(V(:,refAxIdNuc)'*V(:,comp))) *180/pi;
    if ang<minAng
        minAng = ang;
        minAngIndex = comp;
    end
end
refAxIdNuc = minAngIndex;
if minAngIndex == 1
    axId1 = 3;
    axId2 = 2;
elseif minAngIndex == 2
    axId1 = 3;
    axId2 = 1;
else
    axId1 = 2;
    axId2 = 1;
end

end

