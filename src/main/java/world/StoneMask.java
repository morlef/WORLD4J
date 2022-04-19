package world;

import world.common.ForwardRealFFT;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/stonemask.cpp
 */
public record StoneMask(World world) {
    public static void getBaseIndex(double currentPosition, final double[] baseTime, int fs, int[] indexRaw) {
        for (int i = 0; i < baseTime.length; ++i)
            indexRaw[i] = (int) Math.round((currentPosition + baseTime[i]) * fs);
    }

    public static void getMainWindow(double currentPosition, final int[] indexRaw, int baseTimeLength, int fs, double windowLengthInTime, double[] mainWindow) {
        double tmp;
        for (int i = 0; i < baseTimeLength; ++i) {
            tmp = (indexRaw[i] - 1.0) / fs - currentPosition;
            mainWindow[i] = 0.42 +
                    0.5 * Math.cos(2.0 * Math.PI * tmp / windowLengthInTime) +
                    0.08 * Math.cos(4.0 * Math.PI * tmp / windowLengthInTime);
        }
    }

    public static void getDiffWindow(final double[] mainWindow, int baseTimeLength, double[] diffWindow) {
        diffWindow[0] = -mainWindow[1] / 2.0;
        for (int i = 1; i < baseTimeLength - 1; ++i)
            diffWindow[i] = -(mainWindow[i + 1] - mainWindow[i - 1]) / 2.0;
        diffWindow[baseTimeLength - 1] = mainWindow[baseTimeLength - 2] / 2.0;
    }

    public static void getSpectra(World world, final double[] x, int fftSize, final int[] indexRaw, final double[] mainWindow, final double[] diffWindow, int baseTimeLength, final ForwardRealFFT forwardRealFFT, FFT.FFTComplex mainSpectrum, FFT.FFTComplex diffSpectrum) {
        int[] index = new int[baseTimeLength];

        for (int i = 0; i < baseTimeLength; ++i)
            index[i] = Math.max(0, Math.min(x.length - 1, indexRaw[i] - 1));
        for (int i = 0; i < baseTimeLength; ++i)
            forwardRealFFT.getWaveform()[i] = x[index[i]] * mainWindow[i];
        for (int i = baseTimeLength; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;

        world.getFft().execute(forwardRealFFT.getForwardFFT());
        for (int i = 0; i <= fftSize / 2; ++i) {
            mainSpectrum.get()[i][0] = forwardRealFFT.getSpectrum().get()[i][0];
            mainSpectrum.get()[i][1] = forwardRealFFT.getSpectrum().get()[i][1];
        }

        for (int i = 0; i < baseTimeLength; ++i)
            forwardRealFFT.getWaveform()[i] = x[index[i]] * diffWindow[i];
        for (int i = baseTimeLength; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;
        world.getFft().execute(forwardRealFFT.getForwardFFT());
        for (int i = 0; i <= fftSize / 2; ++i) {
            diffSpectrum.get()[i][0] = forwardRealFFT.getSpectrum().get()[i][0];
            diffSpectrum.get()[i][1] = forwardRealFFT.getSpectrum().get()[i][1];
        }
    }

    public static double fixF0(final double[] powerSpectrum, final double[] numeratorI, int fftSize, int fs, double initialF0, int numberOfHarmonics) {
        double[] amplitudeList = new double[numberOfHarmonics];
        double[] instantaneousFrequencyList = new double[numberOfHarmonics];
        int index;
        for (int i = 0; i < numberOfHarmonics; ++i) {
            index = (int) Math.min(Math.round(initialF0 * fftSize / fs * (i + 1)), fftSize / 2);
            instantaneousFrequencyList[i] = powerSpectrum[index] == 0.0 ? 0.0 : (double) (index) * fs / fftSize + numeratorI[index] / powerSpectrum[index] * fs / 2.0 / Math.PI;
            amplitudeList[i] = Math.sqrt(powerSpectrum[index]);
        }
        double denominator = 0.0;
        double numerator = 0.0;
        for (int i = 0; i < numberOfHarmonics; ++i) {
            numerator += amplitudeList[i] * instantaneousFrequencyList[i];
            denominator += amplitudeList[i] * (i + 1);
        }
        return numerator / (denominator + World.MY_SAFE_GUARD_MINIMUM);
    }

    public static double getTentativeF0(final double[] powerSpectrum, final double[] numeratorI, int fftSize, int fs, double initialF0) {
        double tentativeF0 = fixF0(powerSpectrum, numeratorI, fftSize, fs, initialF0, 2);

        if (tentativeF0 <= 0.0 || tentativeF0 > initialF0 * 2) return 0.0;

        return fixF0(powerSpectrum, numeratorI, fftSize, fs, tentativeF0, 6);
    }

    public static double getMeanF0(World world, final double[] x, int fs, double currentPosition, double initialF0, int fftSize, double windowLengthInTime, final double[] baseTime) {
        ForwardRealFFT forwardRealFFT = new ForwardRealFFT(fftSize);
        FFT.FFTComplex mainSpectrum = new FFT.FFTComplex(fftSize);
        FFT.FFTComplex diffSpectrum = new FFT.FFTComplex(fftSize);

        int[] indexRaw = new int[baseTime.length];
        double[] mainWindow = new double[baseTime.length];
        double[] diffWindow = new double[baseTime.length];

        getBaseIndex(currentPosition, baseTime, fs, indexRaw);
        getMainWindow(currentPosition, indexRaw, baseTime.length, fs, windowLengthInTime, mainWindow);
        getDiffWindow(mainWindow, baseTime.length, diffWindow);
        getSpectra(world, x, fftSize, indexRaw, mainWindow, diffWindow, baseTime.length, forwardRealFFT, mainSpectrum, diffSpectrum);

        double[] powerSpectrum = new double[fftSize / 2 + 1];
        double[] numeratorI = new double[fftSize / 2 + 1];
        for (int j = 0; j <= fftSize / 2; ++j) {
            numeratorI[j] = mainSpectrum.get()[j][0] * diffSpectrum.get()[j][1] - mainSpectrum.get()[j][1] * diffSpectrum.get()[j][0];
            powerSpectrum[j] = mainSpectrum.get()[j][0] * mainSpectrum.get()[j][0] + mainSpectrum.get()[j][1] * mainSpectrum.get()[j][1];
        }

        return getTentativeF0(powerSpectrum, numeratorI, fftSize, fs, initialF0);
    }

    public static double getRefinedF0(World world, final double[] x, int fs, double currentPosition, double initialF0) {
        if (initialF0 <= World.FLOOR_F0_STONEMASK || initialF0 > fs / 12.0)
            return 0.0;

        int halfWindowLength = (int) (1.5 * fs / initialF0 + 1.0);
        double windowLengthInTime = (2.0 * halfWindowLength + 1.0) / fs;
        double[] baseTime = new double[halfWindowLength * 2 + 1];
        for (int i = 0; i < halfWindowLength * 2 + 1; i++)
            baseTime[i] = (double) (-halfWindowLength + i) / fs;
        int fftSize = (int) (Math.pow(2.0, 2.0 + (int) (Math.log(halfWindowLength * 2.0 + 1.0) / World.LOG_2)));

        double meanF0 = getMeanF0(world, x, fs, currentPosition, initialF0, fftSize, windowLengthInTime, baseTime);

        if (Math.abs(meanF0 - initialF0) > initialF0 * 0.2) meanF0 = initialF0;

        return meanF0;
    }

    public void stoneMask(final double[] x, int fs, final double[] temporalPositions, final double[] f0, double[] refinedF0) {
        for (int i = 0; i < f0.length; i++)
            refinedF0[i] = getRefinedF0(world, x, fs, temporalPositions[i], f0[i]);
    }
}
