function pb = makePointStruct()
    pb.enterFcn = [];
    pb.exitFcn = [];
    pb.traverseFcn = @(fig, currentPoint) set(fig, 'Pointer', 'hand');
end