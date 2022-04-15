package world.common;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;

@Data
@NoArgsConstructor
@AllArgsConstructor
public class ZeroCrossings {
    double[] negativeIntervalLocations;
    double[] negativeIntervals;
    int numberOfNegatives;
    double[] positiveIntervalLocations;
    double[] positiveIntervals;
    int numberOfPositives;
    double[] peakIntervalLocations;
    double[] peakIntervals;
    int numberOfPeaks;
    double[] dipIntervalLocations;
    double[] dipIntervals;
    int numberOfDips;
}
