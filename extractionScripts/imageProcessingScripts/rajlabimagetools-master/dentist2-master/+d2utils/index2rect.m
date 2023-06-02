function outRect = index2rect(start, stop)
    sz = stop-start;
    outRect = [start(1), start(2), sz(1), sz(2)];
end
