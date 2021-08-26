function [m,cu] = straincurv(midstrain)
for i=1:6
    if (i<4)
    mstrain(i) = midstrain(i);
    else
    curv(i-3)= midstrain(i);
    end   
end
m = mstrain';
cu = curv';
end