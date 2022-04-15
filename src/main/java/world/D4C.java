package world;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.experimental.PackagePrivate;
import world.common.ForwardRealFFT;
import world.util.Utils;

import java.util.Arrays;

@PackagePrivate
public record D4C(World world) {
    public static void setParametersForGetWindowedWaveform(int halfWindowLength, int xLength, double currentPosition, int fs, double currentF0,
                                                    int windowType, double windowLengthRatio, int[] baseIndex, int[] safeIndex, double[] window) {
        for (int i = -halfWindowLength; i <= halfWindowLength; ++i)
            baseIndex[i + halfWindowLength] = i;
        int origin = (int) Math.round(currentPosition * fs + 0.001);
        for (int i = 0; i <= halfWindowLength * 2; ++i)
            safeIndex[i] = Math.min(xLength - 1, Math.max(0, origin + baseIndex[i]));

        double position;
        if (windowType == World.HANNING) {
            for (int i = 0; i <= halfWindowLength * 2; ++i) {
                position = (2.0 * baseIndex[i] / windowLengthRatio) / fs;
                window[i] = 0.5 * Math.cos(Math.PI * position * currentF0) + 0.5;
            }
        } else {
            for (int i = 0; i <= halfWindowLength * 2; ++i) {
                position = (2.0 * baseIndex[i] / windowLengthRatio) / fs;
                window[i] = 0.42 + 0.5 * Math.cos(Math.PI * position * currentF0) +
                        0.08 * Math.cos(Math.PI * position * currentF0 * 2);
            }
        }
    }

    public static void getWindowedWaveform(final double[] x, int fs, double currentF0, double currentPosition, int windowType, double windowLengthRatio, double[] waveform) {
        int halfWindowLength = (int) Math.round(windowLengthRatio * fs / currentF0 / 2.0);

        int[] baseIndex = new int[halfWindowLength * 2 + 1];
        int[] safeIndex = new int[halfWindowLength * 2 + 1];
        double[] window  = new double[halfWindowLength * 2 + 1];

        setParametersForGetWindowedWaveform(halfWindowLength, x.length, currentPosition, fs, currentF0, windowType, windowLengthRatio, baseIndex, safeIndex, window);

        for (int i = 0; i <= halfWindowLength * 2; ++i)
            waveform[i] = x[safeIndex[i]] * window[i] + Math.random() * World.MY_SAFE_GUARD_MINIMUM;

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

    public static void getCentroid(World world, final double[] x, int fs, double currentF0, int fftSize, double currentPosition, final ForwardRealFFT forwardRealFFT, double[] centroid) {
        for (int i = 0; i < fftSize; ++i) forwardRealFFT.getWaveform()[i] = 0.0;
        getWindowedWaveform(x, fs, currentF0, currentPosition, World.BLACKMAN, 4.0, forwardRealFFT.getWaveform());
        double power = 0.0;
        for (int i = 0; i <= Math.round(2.0 * fs / currentF0) * 2; ++i)
            power += forwardRealFFT.getWaveform()[i] * forwardRealFFT.getWaveform()[i];
        for (int i = 0; i <= Math.round(2.0 * fs / currentF0) * 2; ++i)
            forwardRealFFT.getWaveform()[i] /= Math.sqrt(power);

        world.getFft().execute(forwardRealFFT.getForwardFFT());
        double[] tmpReal = new double[fftSize / 2 + 1];
        double[] tmpImag = new double[fftSize / 2 + 1];
        for (int i = 0; i <= fftSize / 2; ++i) {
            tmpReal[i] = forwardRealFFT.getSpectrum().get()[i][0];
            tmpImag[i] = forwardRealFFT.getSpectrum().get()[i][1];
        }

        for (int i = 0; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] *= i + 1.0;
        world.getFft().execute(forwardRealFFT.getForwardFFT());
        for (int i = 0; i <= fftSize / 2; ++i)
            centroid[i] = forwardRealFFT.getSpectrum().get()[i][0] * tmpReal[i] + tmpImag[i] * forwardRealFFT.getSpectrum().get()[i][1];
    }

    public static void getStaticCentroid(World world, final double[] x, int fs, double currentF0, int fftSize, double currentPosition, final ForwardRealFFT forwardRealFFT, double[] staticCentroid) {
        double[] centroid1 = new double[fftSize / 2 + 1];
        double[] centroid2 = new double[fftSize / 2 + 1];

        getCentroid(world, x, fs, currentF0, fftSize, currentPosition - 0.25 / currentF0, forwardRealFFT, centroid1);
        getCentroid(world, x, fs, currentF0, fftSize, currentPosition + 0.25 / currentF0, forwardRealFFT, centroid2);

        for (int i = 0; i <= fftSize / 2; ++i)
            staticCentroid[i] = centroid1[i] + centroid2[i];

        Utils.dcCorrection(staticCentroid, currentF0, fs, fftSize, staticCentroid);
    }

    public static void getSmoothedPowerSpectrum(World world, final double[] x, int fs, double currentF0, int fftSize, double currentPosition, final ForwardRealFFT forwardRealFFT, double[] smoothedPowerSpectrum) {
        for (int i = 0; i < fftSize; ++i) forwardRealFFT.getWaveform()[i] = 0.0;
        getWindowedWaveform(x, fs, currentF0, currentPosition, World.HANNING, 4.0, forwardRealFFT.getWaveform());

        world.getFft().execute(forwardRealFFT.getForwardFFT());
        for (int i = 0; i <= fftSize / 2; ++i)
            smoothedPowerSpectrum[i] = forwardRealFFT.getSpectrum().get()[i][0] * forwardRealFFT.getSpectrum().get()[i][0] + forwardRealFFT.getSpectrum().get()[i][1] * forwardRealFFT.getSpectrum().get()[i][1];
        Utils.dcCorrection(smoothedPowerSpectrum, currentF0, fs, fftSize, smoothedPowerSpectrum);
        Utils.linearSmoothing(smoothedPowerSpectrum, currentF0, fs, fftSize, smoothedPowerSpectrum);
    }

    public static void getStaticGroupDelay(final double[] staticCentroid, final double[] smoothedPowerSpectrum, int fs, double f0, int fftSize, double[] staticGroupDelay) {
        for (int i = 0; i <= fftSize / 2; ++i)
            staticGroupDelay[i] = staticCentroid[i] / smoothedPowerSpectrum[i];
        Utils.linearSmoothing(staticGroupDelay, f0 / 2.0, fs, fftSize, staticGroupDelay);

        double[] smoothedGroupDelay = new double[fftSize / 2 + 1];
        Utils.linearSmoothing(staticGroupDelay, f0, fs, fftSize, smoothedGroupDelay);

        for (int i = 0; i <= fftSize / 2; ++i)
            staticGroupDelay[i] -= smoothedGroupDelay[i];
    }

    public static void getCoarseAperiodicity(World world, final double[] staticGroupDelay, int fs, int fftSize, int numberOfAperiodicities,
                                             final double[] window, final ForwardRealFFT forwardRealFFT, double[] coarseAperiodicity, int s) {
        int boundary = (int) Math.round(fftSize * 8.0 / window.length);
        int halfWindowLength = window.length / 2;

        for (int i = 0; i < fftSize; ++i) forwardRealFFT.getWaveform()[i] = 0.0;

        double[] powerSpectrum = new double[fftSize / 2 + 1];
        int center;
        for (int i = 0; i < numberOfAperiodicities; ++i) {
            center = (int) (World.FREQUENCY_INTERVAL * (i + 1) * fftSize / fs);
            for (int j = 0; j <= halfWindowLength * 2; ++j)
                forwardRealFFT.getWaveform()[j] = staticGroupDelay[center - halfWindowLength + j] * window[j];
            world.getFft().execute(forwardRealFFT.getForwardFFT());
            for (int j = 0 ; j <= fftSize / 2; ++j)
                powerSpectrum[j] = forwardRealFFT.getSpectrum().get()[j][0] * forwardRealFFT.getSpectrum().get()[j][0] + forwardRealFFT.getSpectrum().get()[j][1] * forwardRealFFT.getSpectrum().get()[j][1];
            Arrays.sort(powerSpectrum);
            for (int j = 1 ; j <= fftSize / 2; ++j)
                powerSpectrum[j] += powerSpectrum[j - 1];
            coarseAperiodicity[i + s] = 10 * Math.log10(powerSpectrum[fftSize / 2 - boundary - 1] / powerSpectrum[fftSize / 2]);
        }
    }

    public static double d4cLoveTrainSub(World world, final double[] x, int fs, double currentF0, double currentPosition, int fftSize, int boundary0, int boundary1, int boundary2, ForwardRealFFT forwardRealFFT) {
        double[] powerSpectrum = new double[fftSize];

        int windowLength = (int) (Math.round(1.5 * fs / currentF0) * 2 + 1);
        getWindowedWaveform(x, fs, currentF0, currentPosition, World.BLACKMAN, 3.0, forwardRealFFT.getWaveform());

        for (int i = windowLength; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;
        world.getFft().execute(forwardRealFFT.getForwardFFT());

        for (int i = 0; i <= boundary0; ++i) powerSpectrum[i] = 0.0;
        for (int i = boundary0 + 1; i < fftSize / 2 + 1; ++i)
            powerSpectrum[i] = forwardRealFFT.getSpectrum().get()[i][0] * forwardRealFFT.getSpectrum().get()[i][0] + forwardRealFFT.getSpectrum().get()[i][1] * forwardRealFFT.getSpectrum().get()[i][1];
        for (int i = boundary0; i <= boundary2; ++i)
            powerSpectrum[i] += powerSpectrum[i - 1];

        return powerSpectrum[boundary1] / powerSpectrum[boundary2];
    }

    public static void d4cLoveTrain(World world, final double[] x, int fs, final double[] f0, final double[] temporalPositions, double[] aperiodicity0) {
        double lowestF0 = 40.0;
        int fftSize = (int) (Math.pow(2.0, 1.0 + (int) (Math.log(3.0 * fs / lowestF0 + 1) / World.LOG_2)));
        ForwardRealFFT forwardRealFFT = new ForwardRealFFT();
        Utils.initializeForwardRealFFT(fftSize, forwardRealFFT);

        int boundary0 = (int) (Math.ceil(100.0 * fftSize / fs));
        int boundary1 = (int) (Math.ceil(4000.0 * fftSize / fs));
        int boundary2 = (int) (Math.ceil(7900.0 * fftSize / fs));
        for (int i = 0; i < f0.length; ++i) {
            if (f0[i] == 0.0) {
                aperiodicity0[i] = 0.0;
                continue;
            }
            aperiodicity0[i] = d4cLoveTrainSub(world, x, fs, Math.max(f0[i], lowestF0), temporalPositions[i], fftSize, boundary0, boundary1, boundary2, forwardRealFFT);
        }
    }

    public static void d4cGeneralBody(World world, final double[] x, int fs, double currentF0, int fftSize, double currentPosition,
                                      int numberOfAperiodicities, final double[] window, final ForwardRealFFT forwardRealFFT, double[] coarseAperiodicity, int s) {
        double[] staticCentroid = new double[fftSize / 2 + 1];
        double[] smoothedPowerSpectrum = new double[fftSize / 2 + 1];
        double[] staticGroupDelay = new double[fftSize / 2 + 1];
        getStaticCentroid(world, x, fs, currentF0, fftSize, currentPosition, forwardRealFFT, staticCentroid);
        getSmoothedPowerSpectrum(world, x, fs, currentF0, fftSize, currentPosition, forwardRealFFT, smoothedPowerSpectrum);
        getStaticGroupDelay(staticCentroid, smoothedPowerSpectrum, fs, currentF0, fftSize, staticGroupDelay);

        getCoarseAperiodicity(world, staticGroupDelay, fs, fftSize, numberOfAperiodicities, window, forwardRealFFT, coarseAperiodicity, s);

        for (int i = 0; i < numberOfAperiodicities; ++i)
            coarseAperiodicity[i + s] = Math.min(0.0, coarseAperiodicity[i + s] + (currentF0 - 100) / 50.0);
    }

    public static void initializeAperiodicity(int f0Length, int fftSize, double[][] aperiodicity) {
        for (int i = 0; i < f0Length; ++i)
            for (int j = 0; j < fftSize / 2 + 1; ++j)
                aperiodicity[i][j] = 1.0 - World.MY_SAFE_GUARD_MINIMUM;
    }

    public static void getAperiodicity(final double[] coarseFrequencyAxis, final double[] coarseAperiodicity, final double[] frequencyAxis, double[] aperiodicity) {
        Utils.interp1(coarseFrequencyAxis, coarseAperiodicity, frequencyAxis, aperiodicity);
        for (int i = 0; i <= frequencyAxis.length / 2; ++i)
            aperiodicity[i] = Math.pow(10.0, aperiodicity[i] / 20.0);
    }

    public void d4c(final double[] x, int fs, final double[] temporalPositions, final double[] f0, int fftSize, final D4COption option, double[][] aperiodicity) {

        initializeAperiodicity(f0.length, fftSize, aperiodicity);

        int fftSizeD4C = (int) (Math.pow(2.0, 1.0 + (int) (Math.log(4.0 * fs / World.FLOOR_F0D4C + 1) / World.LOG_2)));

        ForwardRealFFT forwardRealFFT = new ForwardRealFFT();
        Utils.initializeForwardRealFFT(fftSizeD4C, forwardRealFFT);

        int numberOfAperiodicities = (int) (Math.min(World.UPPER_LIMIT, fs / 2.0 - World.FREQUENCY_INTERVAL) / World.FREQUENCY_INTERVAL);
        double[] window =  new double[(int) (World.FREQUENCY_INTERVAL * fftSizeD4C / fs) * 2 + 1];
        Utils.nuttallWindow(window);

        double[] aperiodicity0 = new double[f0.length];
        d4cLoveTrain(world, x, fs, f0, temporalPositions, aperiodicity0);

        double[] coarseAperiodicity = new double[numberOfAperiodicities + 2];
        coarseAperiodicity[0] = -60.0;
        coarseAperiodicity[numberOfAperiodicities + 1] = -World.MY_SAFE_GUARD_MINIMUM;
        double[] coarseFrequencyAxis = new double[numberOfAperiodicities + 2];
        for (int i = 0; i <= numberOfAperiodicities; ++i)
            coarseFrequencyAxis[i] = i * World.FREQUENCY_INTERVAL;
        coarseFrequencyAxis[numberOfAperiodicities + 1] = fs / 2.0;

        double[] frequencyAxis = new double[fftSize / 2 + 1];
        for (int i = 0; i <= fftSize / 2; ++i)
            frequencyAxis[i] = (double) (i) * fs / fftSize;

        for (int i = 0; i < f0.length; ++i) {
            if (f0[i] == 0 || aperiodicity0[i] <= option.threshold) continue;
            d4cGeneralBody(world, x, fs, Math.max(World.FLOOR_F0D4C, f0[i]), fftSizeD4C, temporalPositions[i], numberOfAperiodicities, window, forwardRealFFT, coarseAperiodicity, 1);

            getAperiodicity(coarseFrequencyAxis, coarseAperiodicity, frequencyAxis, aperiodicity[i]);
        }
    }

    public void initializeD4COption(D4COption option) {
        option.threshold = World.THRESHOLD;
    }

    @Data
    @NoArgsConstructor
    @AllArgsConstructor
    public static class D4COption {
        double threshold;
    }
}
