function scanSummary = parseScanSummary(inFileName)
    inFileObj = fopen(inFileName);
    scanSummaryArray = textscan(inFileObj, '%s%q', 'Delimiter', '\t');
    fclose(inFileObj);
    scanSummary = cell2table(scanSummaryArray{2}, 'RowNames', scanSummaryArray{1}');
end