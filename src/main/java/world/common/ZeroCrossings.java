package world.common;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/harvest.cpp
 */
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
