function ret = checkXYPoints( ~, X, Y )
checkedPoints = X.*Y;
checkedPoints(~isfinite(checkedPoints)) = [];
if size(checkedPoints(:),1) > 1
    ret = 1;
else
    ret = 0;
end
end
