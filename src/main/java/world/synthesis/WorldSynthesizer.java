package world.synthesis;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import world.common.ForwardRealFFT;
import world.common.InverseRealFFT;
import world.common.MinimumPhaseAnalysis;

@Data
@AllArgsConstructor
@NoArgsConstructor
public class WorldSynthesizer {
    int fs;
    double framePeriod;
    int bufferSize;
    int numberOfPointers;
    int fftSize;

    double[] buffer;
    int currentPointer;
    int i;

    int[] f0Length;
    int[] f0Origin;
    double[][][] spectrogram;
    double[][][] aperiodicity;

    int currentPointer2;
    int headPointer;
    int synthesizedSample;

    int handoff;
    double handoffPhase;
    double handoffF0;
    int lastLocation;

    int cumulativeFrame;
    int currentFrame;

    double[][] interpolatedVuv;
    double[][] pulseLocations;
    int[][] pulseLocationsIndex;
    int[] numberOfPulses;

    double[] impulseResponse;

    MinimumPhaseAnalysis minimumPhase;
    InverseRealFFT inverseRealFFT;
    ForwardRealFFT forwardRealFFT;
}
