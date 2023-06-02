function size = checkSize(obj)
    tmp = obj;
    S =  whos('tmp');
    size = S.bytes;
    clear tmp
end