function outColors = expressionToColors2(expressionVector)
    outColors = interp1([0 max(expressionVector)], [0.8275 0.8275 0.8275; 0 0 1], expressionVector);
    outColors = single(outColors);
end