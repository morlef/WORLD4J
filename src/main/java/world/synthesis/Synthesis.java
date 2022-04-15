package world.synthesis;

import world.World;
import world.common.ForwardRealFFT;
import world.common.InverseRealFFT;
import world.common.MinimumPhaseAnalysis;
import world.util.Utils;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/synthesis.cpp
 */
public record Synthesis(World world) {
    public static void getNoiseSpectrum(World world, int noiseSize, int fftSize, final ForwardRealFFT forwardRealFFT) {
        double average = 0.0;
        for (int i = 0; i < noiseSize; ++i) {
            forwardRealFFT.getWaveform()[i] = Math.random();
            average += forwardRealFFT.getWaveform()[i];
        }

        average /= noiseSize;
        for (int i = 0; i < noiseSize; ++i)
            forwardRealFFT.getWaveform()[i] -= average;
        for (int i = noiseSize; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;
        world.getFft().execute(forwardRealFFT.getForwardFFT());
    }

    public static void getAperiodicResponse(World world, int noiseSize, int fftSize, final double[] spectrum, final double[] aperiodicRatio, double currentVuv, final ForwardRealFFT forwardRealFFT,
                                            final InverseRealFFT inverseRealFFT, final MinimumPhaseAnalysis minimumPhase, double[] aperiodicResponse) {
        getNoiseSpectrum(world, noiseSize, fftSize, forwardRealFFT);

        if (currentVuv != 0.0)
            for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i)
                minimumPhase.getLogSpectrum()[i] = Math.log(spectrum[i] * aperiodicRatio[i]) / 2.0;
        else
            for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i)
                minimumPhase.getLogSpectrum()[i] = Math.log(spectrum[i]) / 2.0;
        Utils.getMinimumPhaseSpectrum(world, minimumPhase);

        for (int i = 0; i <= fftSize / 2; ++i) {
            inverseRealFFT.getSpectrum().get()[i][0] =
                    minimumPhase.getMinimumPhaseSpectrum().get()[i][0] * forwardRealFFT.getSpectrum().get()[i][0] -
                                    minimumPhase.getMinimumPhaseSpectrum().get()[i][1] * forwardRealFFT.getSpectrum().get()[i][1];
            inverseRealFFT.getSpectrum().get()[i][1] =
                    minimumPhase.getMinimumPhaseSpectrum().get()[i][0] * forwardRealFFT.getSpectrum().get()[i][1] +
                                    minimumPhase.getMinimumPhaseSpectrum().get()[i][1] * forwardRealFFT.getSpectrum().get()[i][0];
        }
        world.getFft().execute(inverseRealFFT.getInverseFFT());
        Utils.fftShift(inverseRealFFT.getWaveform(), aperiodicResponse);
    }

    public static void getSpectrumWithFractionalTimeShift(int fftSize, double coefficient, final InverseRealFFT inverseRealFFT) {
        double re, im, re2, im2;
        for (int i = 0; i <= fftSize / 2; ++i) {
            re = inverseRealFFT.getSpectrum().get()[i][0];
            im = inverseRealFFT.getSpectrum().get()[i][1];
            re2 = Math.cos(coefficient * i);
            im2 = Math.sqrt(1.0 - re2 * re2);

            inverseRealFFT.getSpectrum().get()[i][0] = re * re2 + im * im2;
            inverseRealFFT.getSpectrum().get()[i][1] = im * re2 - re * im2;
        }
    }

    public static void getPeriodicResponse(World world, int fftSize, final double[] spectrum, final double[] aperiodicRatio, double currentVuv, final InverseRealFFT inverseRealFFT,
                                           final MinimumPhaseAnalysis minimumPhase, double fractionalTimeShift, int fs, double[] periodicResponse) {
        if (currentVuv <= 0.5 || aperiodicRatio[0] > 0.999) {
            for (int i = 0; i < fftSize; ++i) periodicResponse[i] = 0.0;
            return;
        }

        for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i)
            minimumPhase.getLogSpectrum()[i] = Math.log(spectrum[i] * (1.0 - aperiodicRatio[i]) + World.MY_SAFE_GUARD_MINIMUM) / 2.0;
        Utils.getMinimumPhaseSpectrum(world, minimumPhase);

        for (int i = 0; i <= fftSize / 2; ++i) {
            inverseRealFFT.getSpectrum().get()[i][0] = minimumPhase.getMinimumPhaseSpectrum().get()[i][0];
            inverseRealFFT.getSpectrum().get()[i][1] = minimumPhase.getMinimumPhaseSpectrum().get()[i][1];
        }

        double coefficient = 2.0 * Math.PI * fractionalTimeShift * fs / fftSize;
        getSpectrumWithFractionalTimeShift(fftSize, coefficient, inverseRealFFT);

        world.getFft().execute(inverseRealFFT.getInverseFFT());
        Utils.fftShift(inverseRealFFT.getWaveform(), periodicResponse);
    }

    public static void getSpectralEnvelope(double currentTime, double framePeriod, int f0Length, final double[][] spectrogram, int fftSize, double[] spectralEnvelope) {
        int currentFrameFloor = Math.min(f0Length - 1, (int) (Math.floor(currentTime / framePeriod)));
        int currentFrameCeil = Math.min(f0Length - 1, (int) (Math.ceil(currentTime / framePeriod)));
        double interpolation = currentTime / framePeriod - currentFrameFloor;

        if (currentFrameFloor == currentFrameCeil)
            for (int i = 0; i <= fftSize / 2; ++i)
                spectralEnvelope[i] = Math.abs(spectrogram[currentFrameFloor][i]);
        else
            for (int i = 0; i <= fftSize / 2; ++i)
                spectralEnvelope[i] = (1.0 - interpolation) * Math.abs(spectrogram[currentFrameFloor][i]) + interpolation * Math.abs(spectrogram[currentFrameCeil][i]);
    }

    public static void getAperiodicRatio(double currentTime, double framePeriod, int f0Length, final double[][] aperiodicity, int fftSize, double[] aperiodicSpectrum) {
        int currentFrameFloor = Math.min(f0Length - 1, (int) (Math.floor(currentTime / framePeriod)));
        int currentFrameCeil = Math.min(f0Length - 1, (int) (Math.ceil(currentTime / framePeriod)));
        double interpolation = currentTime / framePeriod - currentFrameFloor;

        if (currentFrameFloor == currentFrameCeil)
            for (int i = 0; i <= fftSize / 2; ++i)
                aperiodicSpectrum[i] =
                        Math.pow(Utils.getSafeAperiodicity(aperiodicity[currentFrameFloor][i]), 2.0);
        else
            for (int i = 0; i <= fftSize / 2; ++i)
                aperiodicSpectrum[i] = Math.pow((1.0 - interpolation) * Utils.getSafeAperiodicity(aperiodicity[currentFrameFloor][i]) + interpolation * Utils.getSafeAperiodicity(aperiodicity[currentFrameCeil][i]), 2.0);
    }

    public static void getOneFrameSegment(World world, double currentVuv, int noiseSize, final double[][] spectrogram, int fftSize, final double[][] aperiodicity, int f0Length, double framePeriod, double currentTime,
                                          double fractionalTimeShift, int fs, final ForwardRealFFT forwardRealFFT, final InverseRealFFT inverseRealFFT, final MinimumPhaseAnalysis minimumPhase, double[] response) {
        double[] aperiodicResponse = new double[fftSize];
        double[] periodicResponse = new double[fftSize];

        double[] spectralEnvelope = new double[fftSize];
        double[] aperiodicRatio = new double[fftSize];
        getSpectralEnvelope(currentTime, framePeriod, f0Length, spectrogram, fftSize, spectralEnvelope);
        getAperiodicRatio(currentTime, framePeriod, f0Length, aperiodicity, fftSize, aperiodicRatio);

        getPeriodicResponse(world, fftSize, spectralEnvelope, aperiodicRatio, currentVuv, inverseRealFFT, minimumPhase, fractionalTimeShift, fs, periodicResponse);

        getAperiodicResponse(world, noiseSize, fftSize, spectralEnvelope, aperiodicRatio, currentVuv, forwardRealFFT, inverseRealFFT, minimumPhase, aperiodicResponse);

        double sqrtNoiseSize = Math.sqrt(noiseSize);
        for (int i = 0; i < fftSize; ++i)
            response[i] = (periodicResponse[i] * sqrtNoiseSize + aperiodicResponse[i]) / fftSize;
    }

    public static void getTemporalParametersForTimeBase(final double[] f0, int fs, int yLength, double framePeriod, double lowestF0, double[] timeAxis, double[] coarseTimeAxis, double[] coarseF0, double[] coarseVuv) {
        for (int i = 0; i < yLength; ++i)
            timeAxis[i] = i / (double) (fs);
        for (int i = 0; i < f0.length; ++i) {
            coarseTimeAxis[i] = i * framePeriod;
            coarseF0[i] = f0[i] < lowestF0 ? 0.0 : f0[i];
            coarseVuv[i] = coarseF0[i] == 0.0 ? 0.0 : 1.0;
        }
        coarseTimeAxis[f0.length] = f0.length * framePeriod;
        coarseF0[f0.length] = coarseF0[f0.length - 1] * 2 - coarseF0[f0.length - 2];
        coarseVuv[f0.length] = coarseVuv[f0.length - 1] * 2 - coarseVuv[f0.length - 2];
    }

    public static int getPulseLocationsForTimeBase(final double[] interpolatedF0, final double[] timeAxis, int yLength, int fs, double[] pulseLocations, int[] pulseLocationsIndex, double[] pulseLocationsTimeShift) {
        double[] totalPhase = new double[yLength];
        double[] wrapPhase = new double[yLength];
        double[] wrapPhaseAbs = new double[yLength - 1];
        totalPhase[0] = 2.0 * Math.PI * interpolatedF0[0] / fs;
        wrapPhase[0] = (totalPhase[0] % 2.0 * Math.PI);
        for (int i = 1; i < yLength; ++i) {
            totalPhase[i] = totalPhase[i - 1] + 2.0 * Math.PI * interpolatedF0[i] / fs;
            wrapPhase[i] = (totalPhase[i] % 2.0 * Math.PI);
            wrapPhaseAbs[i - 1] = Math.abs(wrapPhase[i] - wrapPhase[i - 1]);
        }

        int numberOfPulses = 0;
        for (int i = 0; i < yLength - 1; ++i) {
            if (wrapPhaseAbs[i] > Math.PI) {
                pulseLocations[numberOfPulses] = timeAxis[i];
                pulseLocationsIndex[numberOfPulses] = i;

                double y1 = wrapPhase[i] - 2.0 * Math.PI;
                double y2 = wrapPhase[i + 1];
                double x = -y1 / (y2 - y1);
                pulseLocationsTimeShift[numberOfPulses] = x / fs;

                ++numberOfPulses;
            }
        }

        return numberOfPulses;
    }

    public static int getTimeBase(final double[] f0, int fs, double framePeriod, int yLength, double lowestF0, double[] pulseLocations, int[] pulseLocationsIndex, double[] pulseLocationsTimeShift, double[] interpolatedVuv) {
        double[] timeAxis = new double[yLength];
        double[] coarseTimeAxis = new double[f0.length + 1];
        double[] coarseF0 = new double[f0.length + 1];
        double[] coarseVuv = new double[f0.length + 1];
        getTemporalParametersForTimeBase(f0, fs, yLength, framePeriod, lowestF0, timeAxis, coarseTimeAxis, coarseF0, coarseVuv);
        double[] interpolatedF0 = new double[yLength];
        Utils.interp1(coarseTimeAxis, coarseF0, timeAxis, interpolatedF0);
        Utils.interp1(coarseTimeAxis, coarseVuv, timeAxis, interpolatedVuv);

        for (int i = 0; i < yLength; ++i) {
            interpolatedVuv[i] = interpolatedVuv[i] > 0.5 ? 1.0 : 0.0;
            interpolatedF0[i] = interpolatedVuv[i] == 0.0 ? World.DEFAULT_F0 : interpolatedF0[i];
        }

        return getPulseLocationsForTimeBase(interpolatedF0, timeAxis, yLength, fs, pulseLocations, pulseLocationsIndex, pulseLocationsTimeShift);
    }

    public void synthesis(final double[] f0, final double[][] spectrogram, final double[][] aperiodicity, int fftSize, double framePeriod, int fs, int y_length, double[] y) {
        double[] impulseResponse = new double[fftSize];

        for (int i = 0; i < y_length; ++i) y[i] = 0.0;

        MinimumPhaseAnalysis minimumPhase = new MinimumPhaseAnalysis();
        Utils.initializeMinimumPhaseAnalysis(fftSize, minimumPhase);
        InverseRealFFT inverseRealFFT = new InverseRealFFT();
        Utils.initializeInverseRealFFT(fftSize, inverseRealFFT);
        ForwardRealFFT forwardRealFFT = new ForwardRealFFT();
        Utils.initializeForwardRealFFT(fftSize, forwardRealFFT);

        double[] pulseLocations = new double[y_length];
        int[] pulseLocationsIndex = new int[y_length];
        double[] pulseLocations_time_shift = new double[y_length];
        double[] interpolatedVuv = new double[y_length];
        int numberOfPulses = getTimeBase(f0, fs, framePeriod / 1000.0, y_length, (double) (fs / fftSize) + 1.0, pulseLocations, pulseLocationsIndex, pulseLocations_time_shift, interpolatedVuv);

        framePeriod /= 1000.0;
        int noiseSize;
        int index, offset, lower_limit, upper_limit;
        for (int i = 0; i < numberOfPulses; ++i) {
            noiseSize = pulseLocationsIndex[Math.min(numberOfPulses - 1, i + 1)] - pulseLocationsIndex[i];

            getOneFrameSegment(world, interpolatedVuv[pulseLocationsIndex[i]], noiseSize, spectrogram, fftSize, aperiodicity, f0.length, framePeriod, pulseLocations[i], pulseLocations_time_shift[i], fs, forwardRealFFT, inverseRealFFT, minimumPhase, impulseResponse);
            offset = pulseLocationsIndex[i] - fftSize / 2 + 1;
            lower_limit = Math.max(0, -offset);
            upper_limit = Math.min(fftSize, y_length - offset);
            for (int j = lower_limit; j < upper_limit; ++j) {
                index = j + offset;
                y[index] += impulseResponse[j];
            }
        }
    }
}
