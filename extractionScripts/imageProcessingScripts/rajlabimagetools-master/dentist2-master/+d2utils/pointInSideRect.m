function l = pointInSideViewRect(rect, point)
    l = rect(1) > point(1,2) && rect(1)+ rect(3) < point(1,2)...
        && rect(2) > point(1,1) && rect(2)+ rect(4) < point(1,1);
end