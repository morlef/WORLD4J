package world.synthesis;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import world.common.ForwardRealFFT;
import world.common.InverseRealFFT;
import world.common.MinimumPhaseAnalysis;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/world/synthrealtime.h
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
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

    public WorldSynthesizer(int fs, double framePeriod, int fftSize, int bufferSize, int numberOfPointers) {
        this.fs = fs;
        this.framePeriod = framePeriod / 1000.0;
        this.bufferSize = bufferSize;
        this.numberOfPointers = numberOfPointers;
        this.fftSize = fftSize;

        this.f0Length = new int[numberOfPointers];
        this.spectrogram = new double[numberOfPointers][][];
        this.aperiodicity = new double[numberOfPointers][][];
        this.interpolatedVuv = new double[numberOfPointers][];
        this.pulseLocations = new double[numberOfPointers][];
        this.pulseLocationsIndex = new int[numberOfPointers][];
        this.numberOfPulses = new int[numberOfPointers];
        this.f0Origin = new int[numberOfPointers];
        for (int i = 0; i < this.numberOfPointers; ++i) {
            this.interpolatedVuv[i] = null;
            this.pulseLocations[i] = null;
            this.pulseLocationsIndex[i] = null;
        }

        this.buffer = new double[bufferSize * 2 + fftSize];
        this.impulseResponse = new double[this.fftSize];

        refreshSynthesizer(this);

        this.minimumPhase = new MinimumPhaseAnalysis(fftSize);
        this.inverseRealFFT = new InverseRealFFT(fftSize);
        this.forwardRealFFT = new ForwardRealFFT(fftSize);
    }

    public void refreshSynthesizer(WorldSynthesizer synth) {
        RealtimeSynthesis.clearRingBuffer(0, synth.numberOfPointers, synth);
        synth.handoffPhase = 0;
        synth.handoffF0 = 0;
        synth.cumulativeFrame = -1;
        synth.lastLocation = 0;

        synth.currentPointer = 0;
        synth.currentPointer2 = 0;
        synth.headPointer = 0;
        synth.handoff = 0;

        synth.i = 0;
        synth.currentFrame = 0;

        synth.synthesizedSample = 0;

        for (int i = 0; i < synth.bufferSize * 2 + synth.fftSize; ++i)
            synth.buffer[i] = 0;
    }
}
