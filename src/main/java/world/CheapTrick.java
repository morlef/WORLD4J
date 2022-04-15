package world;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.experimental.PackagePrivate;
import world.common.ForwardRealFFT;
import world.common.InverseRealFFT;
import world.util.Utils;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/cheaptrick.cpp
 */
@PackagePrivate
public record CheapTrick(World world) {
    public void initialize(final double[] x, int fs, final double[] temporalPositions, final double[] f0, final CheapTrickOption option, double[][] spectrogram) {
        int fftSize = option.fftSize;

        double f0Floor = getF0FloorForCheapTrick(fs, fftSize);
        double[] spectralEnvelope = new double[fftSize];

        ForwardRealFFT forwardRealFFT = new ForwardRealFFT();
        Utils.initializeForwardRealFFT(fftSize, forwardRealFFT);
        InverseRealFFT inverseRealFFT = new InverseRealFFT();
        Utils.initializeInverseRealFFT(fftSize, inverseRealFFT);

        double currentF0;
        for (int i = 0; i < f0.length; ++i) {
            currentF0 = f0[i] <= f0Floor ? World.DEFAULT_F0 : f0[i];
            cheapTrickGeneralBody(x, fs, currentF0, fftSize, temporalPositions[i], option.q1, forwardRealFFT, inverseRealFFT, spectralEnvelope);
            System.arraycopy(spectralEnvelope, 0, spectrogram[i], 0, fftSize / 2 + 1);
        }
    }

    private void initializeCheapTrickOption(int fs, CheapTrickOption option) {
        option.q1 = -0.15;
        option.f0Floor = World.FLOOR_F0;
        option.fftSize = getFFTSizeForCheapTrick(fs, option);
    }

    public int getFFTSizeForCheapTrick(int fs, final CheapTrickOption option) {
        return (int) (Math.pow(2.0, 1.0 + (int) (Math.log(3.0 * fs / option.f0Floor + 1) / World.LOG_2)));
    }

    public double getF0FloorForCheapTrick(int fs, int fftSize) {
        return 3.0 * fs / (fftSize - 3.0);
    }

    public void smoothingWithRecovery(double f0, int fs, int fftSize, double q1, final ForwardRealFFT forwardRealFFT, final InverseRealFFT inverseRealFFT, double[] spectralEnvelope) {
        double[] smoothingLifter = new double[fftSize];
        double[] compensationLifter = new double[fftSize];

        smoothingLifter[0] = 1.0;
        compensationLifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
        double frequency;
        for (int i = 1; i <= forwardRealFFT.getFftSize() / 2; ++i) {
            frequency = (double) (i) / fs;
            smoothingLifter[i] = Math.sin(Math.PI * f0 * frequency) /
                    (Math.PI * f0 * frequency);
            compensationLifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 *
                    Math.cos(2.0 * Math.PI * frequency * f0);
        }

        for (int i = 0; i <= fftSize / 2; ++i)
            forwardRealFFT.getWaveform()[i] = Math.log(forwardRealFFT.getWaveform()[i]);
        for (int i = 1; i < fftSize / 2; ++i)
            forwardRealFFT.getWaveform()[fftSize - i] = forwardRealFFT.getWaveform()[i];
        world.getFft().execute(forwardRealFFT.getForwardFFT());

        for (int i = 0; i <= fftSize / 2; ++i) {
            inverseRealFFT.getSpectrum().get()[i][0] = forwardRealFFT.getSpectrum().get()[i][0] *
                    smoothingLifter[i] * compensationLifter[i] / fftSize;
            inverseRealFFT.getSpectrum().get()[i][1] = 0.0;
        }
        world.getFft().execute(inverseRealFFT.getInverseFFT());

        for (int i = 0; i <= fftSize / 2; ++i)
            spectralEnvelope[i] = Math.exp(inverseRealFFT.getWaveform()[i]);
    }

    public void getPowerSpectrum(int fs, double f0, int fftSize, final ForwardRealFFT forwardRealFFT) {
        int halfWindowLength = (int) Math.round(1.5 * fs / f0);

        for (int i = halfWindowLength * 2 + 1; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;
        world.getFft().execute(forwardRealFFT.getForwardFFT());

        double[] powerSpectrum = forwardRealFFT.getWaveform();
        for (int i = 0; i <= fftSize / 2; ++i)
            powerSpectrum[i] =
                    forwardRealFFT.getSpectrum().get()[i][0] * forwardRealFFT.getSpectrum().get()[i][0] +
                            forwardRealFFT.getSpectrum().get()[i][1] * forwardRealFFT.getSpectrum().get()[i][1];

        Utils.dcCorrection(powerSpectrum, f0, fs, fftSize, powerSpectrum);
    }

    public void setParametersForGetWindowedWaveform(int halfWindowLength, int xLength, double currentPosition, int fs, double currentF0, int[] baseIndex, int[] safeIndex, double[] window) {
        for (int i = -halfWindowLength; i <= halfWindowLength; ++i)
            baseIndex[i + halfWindowLength] = i;
        int origin = (int) Math.round(currentPosition * fs + 0.001);
        for (int i = 0; i <= halfWindowLength * 2; ++i)
            safeIndex[i] =
                    Math.min(xLength - 1, Math.max(0, origin + baseIndex[i]));

        double average = 0.0;
        double position;
        for (int i = 0; i <= halfWindowLength * 2; ++i) {
            position = baseIndex[i] / 1.5 / fs;
            window[i] = 0.5 * Math.cos(Math.PI * position * currentF0) + 0.5;
            average += window[i] * window[i];
        }
        average = Math.sqrt(average);
        for (int i = 0; i <= halfWindowLength * 2; ++i) window[i] /= average;
    }

    public void getWindowedWaveform(final double[] x, int fs, double currentF0, double currentPosition, final ForwardRealFFT forwardRealFFT) {
        int halfWindowLength = (int) Math.round(1.5 * fs / currentF0);

        int[] baseIndex = new int[halfWindowLength * 2 + 1];
        int[] safeIndex = new int[halfWindowLength * 2 + 1];
        double[] window = new double[halfWindowLength * 2 + 1];

        setParametersForGetWindowedWaveform(halfWindowLength, x.length, currentPosition, fs, currentF0, baseIndex, safeIndex, window);

        double[] waveform = forwardRealFFT.getWaveform();
        for (int i = 0; i <= halfWindowLength * 2; ++i)
            waveform[i] = x[safeIndex[i]] * window[i] +
                    Math.random() * World.MY_SAFE_GUARD_MINIMUM;
        double tmpWeight1 = 0;
        double tmpWeight2 = 0;
        for (int i = 0; i <= halfWindowLength * 2; ++i) {
            tmpWeight1 += waveform[i];
            tmpWeight2 += window[i];
        }
        double weightingCoefficient = tmpWeight1 / tmpWeight2;
        for (int i = 0; i <= halfWindowLength * 2; ++i)
            waveform[i] -= window[i] * weightingCoefficient;
    }

    public void addInfinitesimalNoise(final double[] inputSpectrum, int fftSize, double[] outputSpectrum) {
        for (int i = 0; i <= fftSize / 2; ++i)
            outputSpectrum[i] = inputSpectrum[i] + Math.abs(Math.random()) * World.EPS;
    }

    public void cheapTrickGeneralBody(final double[] x, int fs, double currentF0, int fftSize, double currentPosition, double q1,
                                      final ForwardRealFFT forwardRealFFT, final InverseRealFFT inverseRealFFT, double[] spectralEnvelope) {
        getWindowedWaveform(x, fs, currentF0, currentPosition, forwardRealFFT);

        getPowerSpectrum(fs, currentF0, fftSize, forwardRealFFT);

        Utils.linearSmoothing(forwardRealFFT.getWaveform(), currentF0 * 2.0 / 3.0, fs, fftSize, forwardRealFFT.getWaveform());

        addInfinitesimalNoise(forwardRealFFT.getWaveform(), fftSize, forwardRealFFT.getWaveform());

        smoothingWithRecovery(currentF0, fs, fftSize, q1, forwardRealFFT, inverseRealFFT, spectralEnvelope);
    }

    @Data
    @AllArgsConstructor
    @NoArgsConstructor
    public static class CheapTrickOption {
        double q1;
        double f0Floor;
        int fftSize;
    }
}
