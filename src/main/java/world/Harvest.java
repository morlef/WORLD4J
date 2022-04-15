package world;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.experimental.PackagePrivate;
import world.common.ForwardRealFFT;
import world.common.ZeroCrossings;
import world.util.Holder;
import world.util.Utils;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/harvest.cpp
 */
@PackagePrivate
public record Harvest(World world) {
    public static void getWaveformAndSpectrumSub(final double[] x, int decimationRatio, double[] y) {
        if (decimationRatio == 1) {
            System.arraycopy(x, 0, y, 0, x.length);
            return;
        }

        int lag = (int) (Math.ceil(140.0 / decimationRatio) * decimationRatio);
        int newXLength = x.length + lag * 2;
        double[] newY = new double[newXLength];
        for (int i = 0; i < newXLength; ++i) newY[i] = 0.0;
        double[] newX = new double[newXLength];
        for (int i = 0; i < lag; ++i) newX[i] = x[0];
        System.arraycopy(x, 0, newX, lag, lag + x.length - lag);
        for (int i = lag + x.length; i < newXLength; ++i)
            newX[i] = x[x.length - 1];

        Utils.decimate(newX, decimationRatio, newY);
        System.arraycopy(newY, lag / decimationRatio, y, 0, y.length);
    }

    public static void getWaveformAndSpectrum(World world, final double[] x, int fftSize, int decimation_ratio, double[] y, FFT.FFTComplex y_spectrum) {
        for (int i = 0; i < fftSize; ++i) y[i] = 0.0;

        getWaveformAndSpectrumSub(x, decimation_ratio, y);

        double mean_y = 0.0;
        for (double v : y) mean_y += v;
        mean_y /= y.length;
        for (int i = 0; i < y.length; ++i) y[i] -= mean_y;

        FFT.FFTPlan forwardFFT = FFT.r2c(fftSize, y, y_spectrum, FFT.ESTIMATE);
        world.getFft().execute(forwardFFT);
    }

    public static void getFilteredSignal(World world, double boundaryF0, int fftSize, double fs, final FFT.FFTComplex ySpectrum, int yLength, double[] filteredSignal) {
        int filterLengthHalf = (int) Math.round(fs / boundaryF0 * 2.0);
        double[] bandPassFilter = new double[fftSize];
        Utils.nuttallWindow(bandPassFilter);
        for (int i = -filterLengthHalf; i <= filterLengthHalf; ++i)
            bandPassFilter[i + filterLengthHalf] *= Math.cos(2 * Math.PI * boundaryF0 * i / fs);
        for (int i = filterLengthHalf * 2 + 1; i < fftSize; ++i)
            bandPassFilter[i] = 0.0;

        FFT.FFTComplex bandPassFilterSpectrum = new FFT.FFTComplex(fftSize);
        FFT.FFTPlan forwardFFT = FFT.r2c(fftSize, bandPassFilter,
                bandPassFilterSpectrum, FFT.ESTIMATE);
        world.getFft().execute(forwardFFT);

        double tmp = ySpectrum.get()[0][0] * bandPassFilterSpectrum.get()[0][0] - ySpectrum.get()[0][1] * bandPassFilterSpectrum.get()[0][1];
        bandPassFilterSpectrum.get()[0][1] = ySpectrum.get()[0][0] * bandPassFilterSpectrum.get()[0][1] + ySpectrum.get()[0][1] * bandPassFilterSpectrum.get()[0][0];
        bandPassFilterSpectrum.get()[0][0] = tmp;
        for (int i = 1; i <= fftSize / 2; ++i) {
            tmp = ySpectrum.get()[i][0] * bandPassFilterSpectrum.get()[i][0] - ySpectrum.get()[i][1] * bandPassFilterSpectrum.get()[i][1];
            bandPassFilterSpectrum.get()[i][1] = ySpectrum.get()[i][0] * bandPassFilterSpectrum.get()[i][1] + ySpectrum.get()[i][1] * bandPassFilterSpectrum.get()[i][0];
            bandPassFilterSpectrum.get()[i][0] = tmp;
            bandPassFilterSpectrum.get()[fftSize - i - 1][0] = bandPassFilterSpectrum.get()[i][0];
            bandPassFilterSpectrum.get()[fftSize - i - 1][1] = bandPassFilterSpectrum.get()[i][1];
        }

        FFT.FFTPlan inverseFFT = FFT.c2r(fftSize, bandPassFilterSpectrum, filteredSignal, FFT.ESTIMATE);
        world.getFft().execute(inverseFFT);

        int index_bias = filterLengthHalf + 1;
        if (yLength >= 0) System.arraycopy(filteredSignal, index_bias, filteredSignal, 0, yLength);
    }

    public static int checkEvent(int x) {
        return x > 0 ? 1 : 0;
    }

    public static int zeroCrossingEngine(final double[] filteredSignal, int yLength, double fs, double[] intervalLocations, double[] intervals) {
        int[] negativeGoingPoints = new int[yLength];

        for (int i = 0; i < yLength - 1; ++i)
            negativeGoingPoints[i] =
                    0.0 < filteredSignal[i] && filteredSignal[i + 1] <= 0.0 ? i + 1 : 0;
        negativeGoingPoints[yLength - 1] = 0;

        int[] edges = new int[yLength];
        int count = 0;
        for (int i = 0; i < yLength; ++i)
            if (negativeGoingPoints[i] > 0)
                edges[count++] = negativeGoingPoints[i];

        if (count < 2) {
            return 0;
        }

        double[] fineEdges = new double[count];
        for (int i = 0; i < count; ++i)
            fineEdges[i] = edges[i] - filteredSignal[edges[i] - 1] /
                    (filteredSignal[edges[i]] - filteredSignal[edges[i] - 1]);

        for (int i = 0; i < count - 1; ++i) {
            intervals[i] = fs / (fineEdges[i + 1] - fineEdges[i]);
            intervalLocations[i] = (fineEdges[i] + fineEdges[i + 1]) / 2.0 / fs;
        }

        return count - 1;
    }

    public static void getFourZeroCrossingIntervals(double[] filteredSignal, int yLength, double actualFs, ZeroCrossings zeroCrossings) {
        zeroCrossings.setNegativeIntervalLocations(new double[yLength]);
        zeroCrossings.setPositiveIntervalLocations(new double[yLength]);
        zeroCrossings.setPeakIntervalLocations(new double[yLength]);
        zeroCrossings.setDipIntervalLocations(new double[yLength]);
        zeroCrossings.setNegativeIntervals(new double[yLength]);
        zeroCrossings.setPositiveIntervals(new double[yLength]);
        zeroCrossings.setPeakIntervals(new double[yLength]);
        zeroCrossings.setDipIntervals(new double[yLength]);

        zeroCrossings.setNumberOfNegatives(zeroCrossingEngine(filteredSignal, yLength, actualFs, zeroCrossings.getNegativeIntervalLocations(), zeroCrossings.getNegativeIntervals()));

        for (int i = 0; i < yLength; ++i) filteredSignal[i] = -filteredSignal[i];
        zeroCrossings.setNumberOfPositives(zeroCrossingEngine(filteredSignal, yLength, actualFs, zeroCrossings.getPositiveIntervalLocations(), zeroCrossings.getPositiveIntervals()));

        for (int i = 0; i < yLength - 1; ++i) filteredSignal[i] = filteredSignal[i] - filteredSignal[i + 1];
        zeroCrossings.setNumberOfPeaks(zeroCrossingEngine(filteredSignal, yLength - 1, actualFs, zeroCrossings.getPeakIntervalLocations(), zeroCrossings.getPeakIntervals()));

        for (int i = 0; i < yLength - 1; ++i) filteredSignal[i] = -filteredSignal[i];
        zeroCrossings.setNumberOfDips(zeroCrossingEngine(filteredSignal, yLength - 1, actualFs, zeroCrossings.getDipIntervalLocations(), zeroCrossings.getDipIntervals()));
    }

    public static void getF0CandidateContourSub(final double[][] interpolatedF0Set, int f0Length, double f0Floor, double f0Ceil, double boundaryF0, double[] f0Candidate) {
        double upper = boundaryF0 * 1.1;
        double lower = boundaryF0 * 0.9;
        for (int i = 0; i < f0Length; ++i) {
            f0Candidate[i] = (interpolatedF0Set[0][i] +
                    interpolatedF0Set[1][i] + interpolatedF0Set[2][i] +
                    interpolatedF0Set[3][i]) / 4.0;

            if (f0Candidate[i] > upper || f0Candidate[i] < lower ||
                    f0Candidate[i] > f0Ceil || f0Candidate[i] < f0Floor)
                f0Candidate[i] = 0.0;
        }
    }

    public static void getF0CandidateContour(final ZeroCrossings zeroCrossings, double boundaryF0, double f0Floor, double f0Ceil, final double[] temporalPositions, int f0Length, double[] f0Candidate) {
        if (0 == checkEvent(zeroCrossings.getNumberOfNegatives() - 2) *
                checkEvent(zeroCrossings.getNumberOfPositives() - 2) *
                checkEvent(zeroCrossings.getNumberOfPeaks() - 2) *
                checkEvent(zeroCrossings.getNumberOfDips() - 2)) {
            for (int i = 0; i < f0Length; ++i) f0Candidate[i] = 0.0;
            return;
        }

        double[][] interpolatedF0_set = new double[f0Length][4];
        for (int i = 0; i < 4; ++i)
            interpolatedF0_set[i] = new double[f0Length];

        Utils.interp1(zeroCrossings.getNegativeIntervalLocations(), zeroCrossings.getNegativeIntervals(), temporalPositions, interpolatedF0_set[0]);
        Utils.interp1(zeroCrossings.getPositiveIntervalLocations(), zeroCrossings.getPositiveIntervals(), temporalPositions, interpolatedF0_set[1]);
        Utils.interp1(zeroCrossings.getPeakIntervalLocations(), zeroCrossings.getPeakIntervals(), temporalPositions, interpolatedF0_set[2]);
        Utils.interp1(zeroCrossings.getDipIntervalLocations(), zeroCrossings.getDipIntervals(), temporalPositions, interpolatedF0_set[3]);

        getF0CandidateContourSub(interpolatedF0_set, f0Length, f0Floor, f0Ceil, boundaryF0, f0Candidate);
    }

    public static void getF0CandidateFromRawEvent(World world, double boundaryF0, double fs, final FFT.FFTComplex ySpectrum, int yLength, int fftSize,
                                                  double f0Floor, double f0Ceil, final double[] temporalPositions, int f0Length, double[] f0Candidate) {
        double[] filtered_signal = new double[fftSize];
        getFilteredSignal(world, boundaryF0, fftSize, fs, ySpectrum, yLength, filtered_signal);

        ZeroCrossings zeroCrossings = new ZeroCrossings();
        getFourZeroCrossingIntervals(filtered_signal, yLength, fs, zeroCrossings);

        getF0CandidateContour(zeroCrossings, boundaryF0, f0Floor, f0Ceil, temporalPositions, f0Length, f0Candidate);
    }

    public static void getRawF0Candidates(World world, final double[] boundaryF0List, int numberOfBands, double actualFs, int yLength, final double[] temporalPositions, int f0Length,
                                          final FFT.FFTComplex ySpectrum, int fftSize, double f0Floor, double f0Ceil, double[][] rawF0Candidates) {
        for (int i = 0; i < numberOfBands; ++i)
            getF0CandidateFromRawEvent(world, boundaryF0List[i], actualFs, ySpectrum, yLength, fftSize, f0Floor, f0Ceil, temporalPositions, f0Length, rawF0Candidates[i]);
    }

    public static int detectOfficialF0CandidatesSub1(final int[] vuv, int numberOfChannels, int[] st, int[] ed) {
        int numberOfVoicedSections = 0;
        int tmp;
        for (int i = 1; i < numberOfChannels; ++i) {
            tmp = vuv[i] - vuv[i - 1];
            if (tmp == 1) st[numberOfVoicedSections] = i;
            if (tmp == -1) ed[numberOfVoicedSections++] = i;
        }

        return numberOfVoicedSections;
    }

    public static int detectOfficialF0CandidatesSub2(final double[][] rawF0Candidates, int index, int numberOfVoicedSections, final int[] st, final int[] ed, int maxCandidates, double[] f0List) {
        int numberOfCandidates = 0;
        double tmpF0;
        for (int i = 0; i < numberOfVoicedSections; ++i) {
            if (ed[i] - st[i] < 10) continue;

            tmpF0 = 0.0;
            for (int j = st[i]; j < ed[i]; ++j)
                tmpF0 += rawF0Candidates[j][index];
            tmpF0 /= (ed[i] - st[i]);
            f0List[numberOfCandidates++] = tmpF0;
        }

        for (int i = numberOfCandidates; i < maxCandidates; ++i) f0List[i] = 0.0;
        return numberOfCandidates;
    }

    public static int detectOfficialF0Candidates(final double[][] rawF0Candidates, int numberOfChannels, int f0Length, int maxCandidates, double[][] f0Candidates) {
        int numberOfCandidates = 0;

        int[] vuv = new int[numberOfChannels];
        int[] st = new int[numberOfChannels];
        int[] ed = new int[numberOfChannels];
        int numberOfVoicedSections;
        for (int i = 0; i < f0Length; ++i) {
            for (int j = 0; j < numberOfChannels; ++j) vuv[j] = rawF0Candidates[j][i] > 0 ? 1 : 0;
            vuv[0] = vuv[numberOfChannels - 1] = 0;
            numberOfVoicedSections = detectOfficialF0CandidatesSub1(vuv, numberOfChannels, st, ed);
            numberOfCandidates = Math.max(numberOfCandidates, detectOfficialF0CandidatesSub2(rawF0Candidates, i, numberOfVoicedSections, st, ed, maxCandidates, f0Candidates[i]));
        }
        return numberOfCandidates;
    }

    public static void overlapF0Candidates(int f0Length, int numberOfCandidates, double[][] f0Candidates) {
        int n = 3;
        for (int i = 1; i <= n; ++i)
            for (int j = 0; j < numberOfCandidates; ++j) {
                for (int k = i; k < f0Length; ++k) f0Candidates[k][j + (numberOfCandidates * i)] = f0Candidates[k - i][j];
                for (int k = 0; k < f0Length - i; ++k) f0Candidates[k][j + (numberOfCandidates * (i + n))] = f0Candidates[k + i][j];
            }
    }

    public static void getBaseIndex(double currentPosition, final double[] baseTime, int baseTimeLength, double fs, int[] baseIndex) {
        int basic_index = (int) Math.round((currentPosition + baseTime[0]) * fs + 0.001);

        for (int i = 0; i < baseTimeLength; ++i) baseIndex[i] = basic_index + i;
    }

    public static void getMainWindow(double currentPosition, final int[] baseIndex, int baseTimeLength, double fs, double windowLengthInTime, double[] mainWindow) {
        double tmp;
        for (int i = 0; i < baseTimeLength; ++i) {
            tmp = (baseIndex[i] - 1.0) / fs - currentPosition;
            mainWindow[i] = 0.42 + 0.5 * Math.cos(2.0 * Math.PI * tmp / windowLengthInTime) + 0.08 * Math.cos(4.0 * Math.PI * tmp / windowLengthInTime);
        }
    }

    public static void getDiffWindow(final double[] mainWindow, int baseTimeLength, double[] diffWindow) {
        diffWindow[0] = -mainWindow[1] / 2.0;
        for (int i = 1; i < baseTimeLength - 1; ++i) diffWindow[i] = -(mainWindow[i + 1] - mainWindow[i - 1]) / 2.0;
        diffWindow[baseTimeLength - 1] = mainWindow[baseTimeLength - 2] / 2.0;
    }

    public static void getSpectra(World world, final double[] x, int fftSize, final int[] baseIndex, final double[] mainWindow,
                                  final double[] diffWindow, int baseTimeLength, final ForwardRealFFT forwardRealFFT, FFT.FFTComplex mainSpectrum, FFT.FFTComplex diffSpectrum) {
        int[] safe_index = new int[baseTimeLength];

        for (int i = 0; i < baseTimeLength; ++i)
            safe_index[i] = Math.max(0, Math.min(x.length - 1, baseIndex[i] - 1));
        getSpectraSub(world, x, fftSize, mainWindow, baseTimeLength, forwardRealFFT, mainSpectrum, safe_index);

        getSpectraSub(world, x, fftSize, diffWindow, baseTimeLength, forwardRealFFT, diffSpectrum, safe_index);
    }

    private static void getSpectraSub(World world, double[] x, int fftSize, double[] diffWindow, int baseTimeLength, ForwardRealFFT forwardRealFFT, FFT.FFTComplex diffSpectrum, int[] safeIndex) {
        for (int i = 0; i < baseTimeLength; ++i)
            forwardRealFFT.getWaveform()[i] = x[safeIndex[i]] * diffWindow[i];
        for (int i = baseTimeLength; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;

        world.getFft().execute(forwardRealFFT.getForwardFFT());

        for (int i = 0; i <= fftSize / 2; ++i) {
            diffSpectrum.get()[i][0] = forwardRealFFT.getSpectrum().get()[i][0];
            diffSpectrum.get()[i][1] = forwardRealFFT.getSpectrum().get()[i][1];
        }
    }

    public static void fixF0(final double[] powerSpectrum, final double[] numeratorI, int fftSize, double fs, double currentF0, int numberOfHarmonics, Holder<Double> refinedF0, Holder<Double> score) {
        double[] amplitude_list = new double[numberOfHarmonics];
        double[] instantaneous_frequency_list = new double[numberOfHarmonics];

        int index;
        for (int i = 0; i < numberOfHarmonics; ++i) {
            index = (int) Math.round(currentF0 * fftSize / fs * (i + 1));
            instantaneous_frequency_list[i] = powerSpectrum[index] == 0.0 ? 0.0 :
                    (double) (index) * fs / fftSize +
                            numeratorI[index] / powerSpectrum[index] * fs / 2.0 / Math.PI;
            amplitude_list[i] = Math.sqrt(powerSpectrum[index]);
        }
        double denominator = 0.0;
        double numerator = 0.0;

        score.set(0.0);
        for (int i = 0; i < numberOfHarmonics; ++i) {
            numerator += amplitude_list[i] * instantaneous_frequency_list[i];
            denominator += amplitude_list[i] * (i + 1.0);
            score.set(score.get() + Math.abs((instantaneous_frequency_list[i] / (i + 1.0) - currentF0) / currentF0));
        }

        refinedF0.set(numerator / (denominator + World.MY_SAFE_GUARD_MINIMUM));
        score.set(1.0 / ( score.get() / numberOfHarmonics + World.MY_SAFE_GUARD_MINIMUM));
    }

    public static void getMeanF0(World world, final double[] x, double fs, double currentPosition, double currentF0, int fftSize,
                                 double windowLengthInTime, final double[] baseTime, Holder<Double> refinedF0, Holder<Double> refinedScore) {
        ForwardRealFFT forwardRealFFT = new ForwardRealFFT();
        Utils.initializeForwardRealFFT(fftSize, forwardRealFFT);
        FFT.FFTComplex mainSpectrum = new FFT.FFTComplex(fftSize);
        FFT.FFTComplex diffSpectrum = new FFT.FFTComplex(fftSize);

        int[] baseIndex = new int[baseTime.length];
        double[] mainWindow = new double[baseTime.length];
        double[] diffWindow = new double[baseTime.length];

        getBaseIndex(currentPosition, baseTime, baseTime.length, fs, baseIndex);
        getMainWindow(currentPosition, baseIndex, baseTime.length, fs, windowLengthInTime, mainWindow);
        getDiffWindow(mainWindow, baseTime.length, diffWindow);

        getSpectra(world, x, fftSize, baseIndex, mainWindow, diffWindow, baseTime.length, forwardRealFFT, mainSpectrum, diffSpectrum);

        double[] powerSpectrum = new double[fftSize / 2 + 1];
        double[] numeratorI = new double[fftSize / 2 + 1];
        for (int j = 0; j <= fftSize / 2; ++j) {
            numeratorI[j] = mainSpectrum.get()[j][0] * diffSpectrum.get()[j][1] - mainSpectrum.get()[j][1] * diffSpectrum.get()[j][0];
            powerSpectrum[j] = mainSpectrum.get()[j][0] * mainSpectrum.get()[j][0] + mainSpectrum.get()[j][1] * mainSpectrum.get()[j][1];
        }

        int numberOfHarmonics = Math.min((int) (fs / 2.0 / currentF0), 6);
        fixF0(powerSpectrum, numeratorI, fftSize, fs, currentF0, numberOfHarmonics, refinedF0, refinedScore);
    }

    public static void getRefinedF0(World world, final double[] x, double fs, double currentPosition,
                                    double currentF0, double f0Floor, double f0Ceil, Holder<Double> refinedF0, Holder<Double> refinedScore) {
        if (currentF0 <= 0.0) {
            refinedF0.set(0.0);
            refinedScore.set(0.0);
            return;
        }

        int halfWindowLength = (int) (1.5 * fs / currentF0 + 1.0);
        double windowLengthInTime = (2.0 * halfWindowLength + 1.0) / fs;
        double[] baseTime = new double[halfWindowLength * 2 + 1];
        for (int i = 0; i < halfWindowLength * 2 + 1; i++) baseTime[i] = (-halfWindowLength + i) / fs;
        int fftSize = (int) (Math.pow(2.0, 2.0 + (int) (Math.log(halfWindowLength * 2.0 + 1.0) / World.LOG_2)));

        getMeanF0(world, x, fs, currentPosition, currentF0, fftSize, windowLengthInTime, baseTime, refinedF0, refinedScore);

        if (refinedF0.get() < f0Floor || refinedF0.get() > f0Ceil || refinedScore.get() < 2.5){
            refinedF0.set(0.0);
            refinedScore.set(0.0);
        }
    }

    public static void refineF0Candidates(World world, final double[] x, double fs, final double[] temporalPositions, int f0Length,
                                          int maxCandidates, double f0Floor, double f0Ceil, double[][] refinedF0Candidates, double[][] f0Scores) {
        for (int i = 0; i < f0Length; i++)
            for (int j = 0; j < maxCandidates; ++j) {
                Holder<Double> c = new Holder<>(refinedF0Candidates[i][j]);
                Holder<Double> s = new Holder<>(f0Scores[i][j]);

                getRefinedF0(world, x, fs, temporalPositions[i], refinedF0Candidates[i][j], f0Floor, f0Ceil, c, s);

                refinedF0Candidates[i][j] = c.get();
                f0Scores[i][j] = s.get();
            }
    }

    public static double selectBestF0(double referenceF0, final double[] f0Candidates, int numberOfCandidates, double allowedRange, Holder<Double> bestError) {
        double bestF0 = 0.0;
        bestError.set(allowedRange);

        double tmp;
        for (int i = 0; i < numberOfCandidates; ++i) {
            tmp = Math.abs(referenceF0 - f0Candidates[i]) / referenceF0;
            if (tmp > bestError.get()) continue;
            bestF0 = f0Candidates[i];
            bestError.set(tmp);
        }

        return bestF0;
    }

    public static void removeUnreliableCandidatesSub(int i, int j, final double[][] tmpF0Candidates, int numberOfCandidates, double[][] f0Candidates, double[][] f0Scores) {
        double referenceF0 = f0Candidates[i][j];
        Holder<Double> error1 = new Holder<>(0.0), error2 = new Holder<>(0.0);
        double minError;
        double threshold = 0.05;
        if (referenceF0 == 0) return;
        selectBestF0(referenceF0, tmpF0Candidates[i + 1], numberOfCandidates, 1.0, error1);
        selectBestF0(referenceF0, tmpF0Candidates[i - 1], numberOfCandidates, 1.0, error2);
        minError = Math.min(error1.get(), error2.get());
        if (minError <= threshold) return;
        f0Candidates[i][j] = 0;
        f0Scores[i][j] = 0;
    }

    public static void removeUnreliableCandidates(int f0Length, int numberOfCandidates, double[][] f0Candidates, double[][] f0Scores) {
        double[][] tmpF0Candidates = new double[numberOfCandidates][f0Length];
        for (int i = 0; i < f0Length; ++i)
            tmpF0Candidates[i] = new double[numberOfCandidates];
        for (int i = 0; i < f0Length; ++i)
            System.arraycopy(f0Candidates[i], 0, tmpF0Candidates[i], 0, numberOfCandidates);

        for (int i = 1; i < f0Length - 1; ++i)
            for (int j = 0; j < numberOfCandidates; ++j)
                removeUnreliableCandidatesSub(i, j, tmpF0Candidates, numberOfCandidates, f0Candidates, f0Scores);
    }

    public static void searchF0Base(final double[][] f0Candidates, final double[][] f0Scores, int f0Length, int numberOfCandidates, double[] baseF0Contour) {
        double tmpBestScore;
        for (int i = 0; i < f0Length; ++i) {
            baseF0Contour[i] = tmpBestScore = 0.0;
            for (int j = 0; j < numberOfCandidates; ++j)
                if (f0Scores[i][j] > tmpBestScore) {
                    baseF0Contour[i] = f0Candidates[i][j];
                    tmpBestScore = f0Scores[i][j];
                }
        }
    }

    public static void fixStep1(final double[] f0Base, double allowedRange, double[] f0Step1) {
        for (int i = 0; i < f0Base.length; ++i) f0Step1[i] = 0.0;
        double referenceF0;
        for (int i = 2; i < f0Base.length; ++i) {
            if (f0Base[i] == 0.0) continue;
            referenceF0 = f0Base[i - 1] * 2 - f0Base[i - 2];
            f0Step1[i] = Math.abs((f0Base[i] - referenceF0) / referenceF0) > allowedRange && Math.abs((f0Base[i] - f0Base[i - 1])) / f0Base[i - 1] > allowedRange ? 0.0 : f0Base[i];
        }
    }

    public static int getBoundaryList(final double[] f0, int[] boundaryList) {
        int numberOfBoundaries = 0;
        int[] vuv = new int[f0.length];
        for (int i = 0; i < f0.length; ++i)
            vuv[i] = f0[i] > 0 ? 1 : 0;
        vuv[0] = vuv[f0.length - 1] = 0;

        for (int i = 1; i < f0.length; ++i)
            if (vuv[i] - vuv[i - 1] != 0) {
                boundaryList[numberOfBoundaries] = i - numberOfBoundaries % 2;
                numberOfBoundaries++;
            }
        return numberOfBoundaries;
    }

    public static void fixStep2(final double[] f0Step1, int voiceRangeMinimum, double[] f0Step2) {
        System.arraycopy(f0Step1, 0, f0Step2, 0, f0Step1.length);
        int[] boundaryList = new int[f0Step1.length];
        int numberOfBoundaries = getBoundaryList(f0Step1, boundaryList);

        for (int i = 0; i < numberOfBoundaries / 2; ++i) {
            if (boundaryList[i * 2 + 1] - boundaryList[i * 2] >= voiceRangeMinimum) continue;
            for (int j = boundaryList[i * 2]; j <= boundaryList[(i * 2) + 1]; ++j) f0Step2[j] = 0.0;
        }
    }

    public static void getMultiChannelF0(final double[] f0, final int[] boundaryList, int numberOfBoundaries, double[][] multiChannelF0) {
        for (int i = 0; i < numberOfBoundaries / 2; ++i) {
            for (int j = 0; j < boundaryList[i * 2]; ++j)
                multiChannelF0[i][j] = 0.0;
            if (boundaryList[i * 2 + 1] + 1 - boundaryList[i * 2] >= 0)
                System.arraycopy(f0, boundaryList[i * 2], multiChannelF0[i], boundaryList[i * 2], boundaryList[i * 2 + 1] + 1 - boundaryList[i * 2]);
            for (int j = boundaryList[i * 2 + 1] + 1; j < f0.length; ++j)
                multiChannelF0[i][j] = 0.0;
        }
    }

    public static int extendF0(int origin, int lastPoint, int shift, final double[][] f0Candidates, int numberOfCandidates, double allowedRange, double[] extendedF0) {
        int threshold = 4;
        double tmpF0 = extendedF0[origin];
        int shiftedOrigin = origin;

        int distance = Math.abs(lastPoint - origin);
        int[] indexList = new int[distance + 1];
        for (int i = 0; i <= distance; ++i) indexList[i] = origin + shift * i;

        int count = 0;
        double dummy = 0.0;
        for (int i = 0; i <= distance; ++i) {
            extendedF0[indexList[i] + shift] = selectBestF0(tmpF0, f0Candidates[indexList[i] + shift], numberOfCandidates, allowedRange, new Holder<>(dummy));
            if (extendedF0[indexList[i] + shift] == 0.0) {
                count++;
            } else {
                tmpF0 = extendedF0[indexList[i] + shift];
                count = 0;
                shiftedOrigin = indexList[i] + shift;
            }
            if (count == threshold) break;
        }
        return shiftedOrigin;
    }

    public static void swap(int index1, int index2, double[][] f0, int[] boundary) {
        double[] tmpPointer;
        int tmpIndex;
        tmpPointer = f0[index1];
        f0[index1] = f0[index2];
        f0[index2] = tmpPointer;
        tmpIndex = boundary[index1 * 2];
        boundary[index1 * 2] = boundary[index2 * 2];
        boundary[index2 * 2] = tmpIndex;
        tmpIndex = boundary[index1 * 2 + 1];
        boundary[index1 * 2 + 1] = boundary[index2 * 2 + 1];
        boundary[index2 * 2 + 1] = tmpIndex;
    }

    public static int extendSub(final double[][] extendedF0, final int[] boundaryList, int numberOfSections, double[][] selectedExtendedF0, int[] selectedBoundaryList) {
        double threshold = 2200.0;
        int count = 0;
        double meanF0 = 0.0;
        int st, ed;
        for (int i = 0; i < numberOfSections; ++i) {
            st = boundaryList[i * 2];
            ed = boundaryList[i * 2 + 1];
            for (int j = st; j < ed; ++j) meanF0 += extendedF0[i][j];
            meanF0 /= ed - st;
            if (threshold / meanF0 < ed - st)
                swap(count++, i, selectedExtendedF0, selectedBoundaryList);
        }
        return count;
    }

    public static int extend(final double[][] multi_channel_f0, int number_of_sections, int f0Length, final int[] boundary_list,
                             final double[][] f0_candidates, int number_of_candidates, double allowed_range, double[][] extended_f0, int[] shifted_boundary_list) {
        int threshold = 100;
        for (int i = 0; i < number_of_sections; ++i) {
            shifted_boundary_list[i * 2 + 1] = extendF0(boundary_list[i * 2 + 1], Math.min(f0Length - 2, boundary_list[i * 2 + 1] + threshold), 1, f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
            shifted_boundary_list[i * 2] = extendF0(boundary_list[i * 2], Math.max(1, boundary_list[i * 2] - threshold), -1, f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
        }

        return extendSub(multi_channel_f0, shifted_boundary_list, number_of_sections, extended_f0, shifted_boundary_list);
    }

    public static void makeSortedOrder(final int[] boundaryList, int numberOfSections, int[] order) {
        for (int i = 0; i < numberOfSections; ++i) order[i] = i;
        int tmp;
        for (int i = 1; i < numberOfSections; ++i)
            for (int j = i - 1; j >= 0; --j)
                if (boundaryList[order[j] * 2] > boundaryList[order[i] * 2]) {
                    tmp = order[i];
                    order[i] = order[j];
                    order[j] = tmp;
                } else {
                    break;
                }
    }

    public static double searchScore(double f0, final double[] f0Candidates, final double[] f0Scores, int numberOfCandidates) {
        double score = 0.0;
        for (int i = 0; i < numberOfCandidates; ++i)
            if (f0 == f0Candidates[i] && score < f0Scores[i]) score = f0Scores[i];
        return score;
    }

    public static int mergeF0Sub(final double[] f0_1, int st1, int ed1, final double[] f0_2, int st2, int ed2, final double[][] f0Candidates, final double[][] f0Scores, int numberOfCandidates, double[] mergedF0) {
        if (st1 <= st2 && ed1 >= ed2) return ed1;

        double score1 = 0.0;
        double score2 = 0.0;
        for (int i = st2; i <= ed1; ++i) {
            score1 += searchScore(f0_1[i], f0Candidates[i], f0Scores[i], numberOfCandidates);
            score2 += searchScore(f0_2[i], f0Candidates[i], f0Scores[i], numberOfCandidates);
        }
        if (score1 > score2)
            if (ed2 + 1 - ed1 >= 0) System.arraycopy(f0_2, ed1, mergedF0, ed1, ed2 + 1 - ed1);
        else if (ed2 + 1 - st2 >= 0) System.arraycopy(f0_2, st2, mergedF0, st2, ed2 + 1 - st2);

        return ed2;
    }

    public static void mergeF0(final double[][] multiChannelF0, int[] boundaryList, int numberOfChannels, int f0Length,
                               final double[][] f0Candidates, final double[][] f0Scores, int numberOfCandidates, double[] mergedF0) {
        int[] order = new int[numberOfChannels];
        makeSortedOrder(boundaryList, numberOfChannels, order);

        if (f0Length >= 0) System.arraycopy(multiChannelF0[0], 0, mergedF0, 0, f0Length);

        for (int i = 1; i < numberOfChannels; ++i)
            if (boundaryList[order[i] * 2] - boundaryList[1] > 0) {
                if (boundaryList[order[i] * 2 + 1] + 1 - boundaryList[order[i] * 2] >= 0)
                    System.arraycopy(multiChannelF0[order[i]], boundaryList[order[i] * 2], mergedF0, boundaryList[order[i] * 2], boundaryList[order[i] * 2 + 1] + 1 - boundaryList[order[i] * 2]);
                boundaryList[0] = boundaryList[order[i] * 2];
                boundaryList[1] = boundaryList[order[i] * 2 + 1];
            } else {
                boundaryList[1] = mergeF0Sub(mergedF0, boundaryList[0], boundaryList[1], multiChannelF0[order[i]], boundaryList[order[i] * 2],
                        boundaryList[order[i] * 2 + 1], f0Candidates, f0Scores, numberOfCandidates, mergedF0);
            }
    }

    public static void fixStep3(final double[] f0Step2, int numberOfCandidates, final double[][] f0Candidates, double allowedRange, final double[][] f0Scores, double[] f0Step3) {
        System.arraycopy(f0Step2, 0, f0Step3, 0, f0Step2.length);
        int[] boundaryList = new int[f0Step2.length];
        int numberOfBoundaries = getBoundaryList(f0Step2, boundaryList);

        double[][] multiChannelF0 = new double[f0Step2.length][numberOfBoundaries / 2];
        for (int i = 0; i < numberOfBoundaries / 2; ++i)
            multiChannelF0[i] = new double[f0Step2.length];
        getMultiChannelF0(f0Step2, boundaryList, numberOfBoundaries, multiChannelF0);

        int numberOfChannels = extend(multiChannelF0, numberOfBoundaries / 2, f0Step2.length, boundaryList, f0Candidates, numberOfCandidates, allowedRange, multiChannelF0, boundaryList);

        if (numberOfChannels != 0)
            mergeF0(multiChannelF0, boundaryList, numberOfChannels, f0Step2.length, f0Candidates, f0Scores, numberOfCandidates, f0Step3);
    }

    public static void fixStep4(final double[] f0Step3, int threshold, double[] f0Step4) {
        System.arraycopy(f0Step3, 0, f0Step4, 0, f0Step3.length);
        int[] boundaryList = new int[f0Step3.length];
        int numberOfBoundaries = getBoundaryList(f0Step3, boundaryList);

        int distance;
        double tmp0, tmp1, coefficient;
        int count;
        for (int i = 0; i < numberOfBoundaries / 2 - 1; ++i) {
            distance = boundaryList[(i + 1) * 2] - boundaryList[i * 2 + 1] - 1;
            if (distance >= threshold) continue;
            tmp0 = f0Step3[boundaryList[i * 2 + 1]] + 1;
            tmp1 = f0Step3[boundaryList[(i + 1) * 2]] - 1;
            coefficient = (tmp1 - tmp0) / (distance + 1.0);
            count = 1;
            for (int j = boundaryList[i * 2 + 1] + 1;
                 j <= boundaryList[(i + 1) * 2] - 1; ++j)
                f0Step4[j] = tmp0 + coefficient * count++;
        }
    }

    public static void fixF0Contour(final double[][] f0_candidates, final double[][] f0_scores, int f0Length, int number_of_candidates, double[] best_f0_contour) {
        double[] tmp_f0_contour1 = new double[f0Length];
        double[] tmp_f0_contour2 = new double[f0Length];

        searchF0Base(f0_candidates, f0_scores, f0Length,
                number_of_candidates, tmp_f0_contour1);
        fixStep1(tmp_f0_contour1, 0.008, tmp_f0_contour2);
        fixStep2(tmp_f0_contour2, 6, tmp_f0_contour1);
        fixStep3(tmp_f0_contour1, number_of_candidates, f0_candidates, 0.18, f0_scores, tmp_f0_contour2);
        fixStep4(tmp_f0_contour2, 9, best_f0_contour);
    }

    public static void filteringF0(final double[] a, final double[] b, double[] x, int st, int ed, double[] y) {
        double[] w ={0.0, 0.0};
        double wt;
        double[] tmpX = new double[x.length];

        for (int i = 0; i < st; ++i) x[i] = x[st];
        for (int i = ed + 1; i < x.length; ++i) x[i] = x[ed];

        for (int i = 0; i < x.length; ++i) {
            wt = x[i] + a[0] * w[0] + a[1] * w[1];
            tmpX[x.length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
            w[1] = w[0];
            w[0] = wt;
        }

        w[0] = w[1] = 0.0;
        for (int i = 0; i < x.length; ++i) {
            wt = tmpX[i] + a[0] * w[0] + a[1] * w[1];
            y[x.length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
            w[1] = w[0];
            w[0] = wt;
        }
    }

    public static void smoothF0Contour(final double[] f0, double[] smoothedF0) {
        final double[] b = {0.0078202080334971724, 0.015640416066994345};
        final double[] a = {1.7347257688092754, -0.76600660094326412};
        int lag = 300;
        int newF0Length = f0.length + lag * 2;
        double[] f0Contour = new double[newF0Length];
        for (int i = 0; i < lag; ++i) f0Contour[i] = 0.0;
        System.arraycopy(f0, 0, f0Contour, lag, lag + f0.length - lag);
        for (int i = lag + f0.length; i < newF0Length; ++i) f0Contour[i] = 0.0;

        int[] boundaryList = new int[newF0Length];
        int numberOfBoundaries = getBoundaryList(f0Contour, boundaryList);
        double[][] multiChannelF0 = new double[newF0Length][numberOfBoundaries / 2];
        for (int i = 0; i < numberOfBoundaries / 2; ++i)
            multiChannelF0[i] = new double[newF0Length];
        getMultiChannelF0(f0Contour,  boundaryList, numberOfBoundaries, multiChannelF0);

        for (int i = 0; i < numberOfBoundaries / 2; ++i) {
            filteringF0(a, b, multiChannelF0[i], boundaryList[i * 2], boundaryList[i * 2 + 1], f0Contour);
            if (boundaryList[i * 2 + 1] + 1 - boundaryList[i * 2] >= 0)
                System.arraycopy(f0Contour, boundaryList[i * 2], smoothedF0, boundaryList[i * 2] - lag, boundaryList[i * 2 + 1] + 1 - boundaryList[i * 2]);
        }
    }

    public static int harvestGeneralBodySub(World world, final double[] boundaryF0List, int numberOfChannels, int f0Length, double actualFs, int yLength,
                                            final double[] temporalPositions, final FFT.FFTComplex ySpectrum, int fftSize, double f0Floor, double f0Ceil, int maxCandidates, double[][] f0Candidates) {
        double[][] rawF0Candidates = new double[f0Length][numberOfChannels];
        for (int i = 0; i < numberOfChannels; ++i)
            rawF0Candidates[i] = new double[f0Length];

        getRawF0Candidates(world, boundaryF0List, numberOfChannels, actualFs, yLength, temporalPositions, f0Length, ySpectrum, fftSize, f0Floor, f0Ceil, rawF0Candidates);

        int number_of_candidates = detectOfficialF0Candidates(rawF0Candidates, numberOfChannels, f0Length, maxCandidates, f0Candidates);

        overlapF0Candidates(f0Length, number_of_candidates, f0Candidates);

        return number_of_candidates;
    }

    public static void harvestGeneralBody(World world, final double[] x, int fs, int framePeriod, double f0Floor, double f0Ceil,
                                          double channelsInOctave, int speed, double[] temporalPositions, double[] f0) {
        double adjustedF0Floor = f0Floor * 0.9;
        double adjustedF0Ceil = f0Ceil * 1.1;
        int numberOfChannels = (int) (Math.log(adjustedF0Ceil / adjustedF0Floor) / World.LOG_2 * channelsInOctave);
        double[] boundaryF0List = new double[numberOfChannels];
        for (int i = 0; i < numberOfChannels; ++i)
            boundaryF0List[i] = adjustedF0Floor * Math.pow(2.0, (i + 1) / channelsInOctave);

        int decimationRatio = Math.max(Math.min(speed, 12), 1);
        int yLength = (int) (Math.ceil((double) (x.length) / decimationRatio));
        double actualFs = (double) (fs) / decimationRatio;
        int fftSize = Utils.getSuitableFFTSize(yLength + 5 + 2 * (int) (2.0 * actualFs / boundaryF0List[0]));

        double[] y = new double[fftSize];
        FFT.FFTComplex ySpectrum = new FFT.FFTComplex(fftSize);
        getWaveformAndSpectrum(world, x, fftSize, decimationRatio, y, ySpectrum);

        int f0Length = world.getHarvest().getSamplesForHarvest(fs, x.length, framePeriod);
        for (int i = 0; i < f0Length; ++i) {
            temporalPositions[i] = i * framePeriod / 1000.0;
            f0[i] = 0.0;
        }

        int overlapParameter = 7;
        int maxCandidates = (int) (Math.round(numberOfChannels / 10.0) * overlapParameter);
        double[][] f0Candidates = new double[maxCandidates][f0Length];
        double[][] f0CandidatesScore = new double[maxCandidates][f0Length];
        for (int i = 0; i < f0Length; ++i) {
            f0Candidates[i] = new double[maxCandidates];
            f0CandidatesScore[i] = new double[maxCandidates];
        }

        int numberOfCandidates = harvestGeneralBodySub(world, boundaryF0List, numberOfChannels, f0Length, actualFs, yLength,
                temporalPositions, ySpectrum, fftSize, f0Floor, f0Ceil, maxCandidates, f0Candidates) * overlapParameter;

        refineF0Candidates(world, y, actualFs, temporalPositions, f0Length,
                numberOfCandidates, f0Floor, f0Ceil, f0Candidates,
                f0CandidatesScore);
        removeUnreliableCandidates(f0Length, numberOfCandidates,
                f0Candidates, f0CandidatesScore);

        double[] bestF0Contour = new double[f0Length];
        fixF0Contour(f0Candidates, f0CandidatesScore, f0Length,
                numberOfCandidates, bestF0Contour);
        smoothF0Contour(bestF0Contour, f0);
    }

    public int getSamplesForHarvest(int fs, int xLength, double framePeriod) {
        return (int) (1000.0 * xLength / fs / framePeriod) + 1;
    }

    public void harvest(final double[] x, int fs, final HarvestOption option, double[] temporalPositions, double[] f0) {
        double targetFs = 8000.0;
        int dimensionRatio = (int) Math.round(fs / targetFs);
        double channelsInOctave = 40;

        if (option.framePeriod == 1.0) {
            harvestGeneralBody(world, x, fs, 1, option.f0Floor, option.f0Ceil, channelsInOctave, dimensionRatio, temporalPositions, f0);
            return;
        }

        int basicFramePeriod = 1;
        int basicF0Length = getSamplesForHarvest(fs, x.length, basicFramePeriod);
        double[] basicF0 = new double[basicF0Length];
        double[] basicTemporalPositions = new double[basicF0Length];
        harvestGeneralBody(world, x, fs, basicFramePeriod, option.f0Floor, option.f0Ceil, channelsInOctave, dimensionRatio, basicTemporalPositions, basicF0);

        int f0Length = getSamplesForHarvest(fs, x.length, option.framePeriod);
        for (int i = 0; i < f0Length; ++i) {
            temporalPositions[i] = i * option.framePeriod / 1000.0;
            f0[i] = basicF0[(int) Math.min(basicF0Length - 1, Math.round(temporalPositions[i] * 1000.0))];
        }
    }

    public void initializeHarvestOption(HarvestOption option) {
        option.f0Ceil = World.CEIL_F0;
        option.f0Floor = World.FLOOR_F0;
        option.framePeriod = 5;
    }

    @Data
    @NoArgsConstructor
    @AllArgsConstructor
    public static class HarvestOption {
        double f0Floor;
        double f0Ceil;
        double framePeriod;
    }
}
