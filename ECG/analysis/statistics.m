% Extract the reference data
refBPM = resultsTable2.BPM(1);
refStdBPM = resultsTable2.StdPeakInterval(1);
refPeakHeight = resultsTable2.AveragePeakHeight(1);
refStdPeakHeight = resultsTable2.StdPeakHeight(1);
refTotalPeaks = resultsTable2.TotalPeaks(1);

for i = 2:height(resultsTable2)
    % Extract current row data
    currBPM = resultsTable2.BPM(i);
    currStdBPM = resultsTable2.StdPeakInterval(i);
    currPeakHeight = resultsTable2.AveragePeakHeight(i);
    currStdPeakHeight = resultsTable2.StdPeakHeight(i);
    currTotalPeaks = resultsTable2.TotalPeaks(i);
    
    % Calculate t-test for BPM
    pooledStdBPM = sqrt(((refStdBPM^2 * (refTotalPeaks - 1)) + (currStdBPM^2 * (currTotalPeaks - 1))) / (refTotalPeaks + currTotalPeaks - 2));
    t_bpm = (refBPM - currBPM) / (pooledStdBPM * sqrt(1/refTotalPeaks + 1/currTotalPeaks));
    df_bpm = refTotalPeaks + currTotalPeaks - 2;
    p_bpm = 2 * tcdf(-abs(t_bpm), df_bpm);
    
    % Calculate t-test for AveragePeakHeight
    pooledStdPeakHeight = sqrt(((refStdPeakHeight^2 * (refTotalPeaks - 1)) + (currStdPeakHeight^2 * (currTotalPeaks - 1))) / (refTotalPeaks + currTotalPeaks - 2));
    t_peak = (refPeakHeight - currPeakHeight) / (pooledStdPeakHeight * sqrt(1/refTotalPeaks + 1/currTotalPeaks));
    df_peak = refTotalPeaks + currTotalPeaks - 2;
    p_peak = 2 * tcdf(-abs(t_peak), df_peak);
    
    % Display results
    fprintf('Comparison between 1000 Hz and %d Hz:\n', resultsTable2.Frequency(i));
    fprintf('BPM: p-value = %.4e\n', p_bpm);
    fprintf('Average Peak Height: p-value = %.4e\n\n', p_peak);
end
