function l = pointInSideAxes(limits, point) %limits = [xmin xmax ymin ymax]
    l = point(1,1) > limits(1) && point(1,1) < limits(2)...
        && point(1,2) > limits(3) && point(1,2) < limits(4);
end