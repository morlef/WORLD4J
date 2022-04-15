package world;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.experimental.PackagePrivate;
import world.common.ZeroCrossings;
import world.util.Holder;
import world.util.Utils;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/dio.cpp
 */
@PackagePrivate
public record Dio(World world) {
    public static void designLowCutFilter(int n, int fftSize, double[] lowCutFilter) {
        for (int i = 1; i <= n; ++i)
            lowCutFilter[i - 1] = 0.5 - 0.5 * Math.cos(i * 2.0 * Math.PI / (n + 1));
        for (int i = n; i < fftSize; ++i)
            lowCutFilter[i] = 0.0;
        double sumOfAmplitude = 0.0;
        for (int i = 0; i < n; ++i)
            sumOfAmplitude += lowCutFilter[i];
        for (int i = 0; i < n; ++i)
            lowCutFilter[i] = -lowCutFilter[i] / sumOfAmplitude;
        if ((n - 1) / 2 >= 0) System.arraycopy(lowCutFilter, 0, lowCutFilter, fftSize - (n - 1) / 2, (n - 1) / 2);
        if (n >= 0) System.arraycopy(lowCutFilter, (n - 1) / 2, lowCutFilter, 0, n);
        lowCutFilter[0] += 1.0;
    }

    public static void getSpectrumForEstimation(World world, final double[] x, int yLength, double actualFs, int fftSize, int decimationRatio, FFT.FFTComplex ySpectrum) {
        double[] y = new double[fftSize];

        for (int i = 0; i < fftSize; ++i)
            y[i] = 0.0;

        if (decimationRatio != 1)
            Utils.decimate(x, decimationRatio, y);
        else
            System.arraycopy(x, 0, y, 0, x.length);

        double meanY = 0.0;
        for (int i = 0; i < yLength; ++i)
            meanY += y[i];
        meanY /= yLength;
        for (int i = 0; i < yLength; ++i)
            y[i] -= meanY;
        for (int i = yLength; i < fftSize; ++i)
            y[i] = 0.0;

        FFT.FFTPlan forwardFFT = FFT.r2c(fftSize, y, ySpectrum, FFT.ESTIMATE);
        world.getFft().execute(forwardFFT);

        int cutoffInSample = (int) Math.round(actualFs / World.CUT_OFF);
        designLowCutFilter(cutoffInSample * 2 + 1, fftSize, y);

        FFT.FFTComplex filterSpectrum = new FFT.FFTComplex(fftSize);
        forwardFFT.cOut = filterSpectrum;
        world.getFft().execute(forwardFFT);

        double tmp;
        for (int i = 0; i <= fftSize / 2; ++i) {
            tmp = ySpectrum.get()[i][0] * filterSpectrum.get()[i][0] - ySpectrum.get()[i][1] * filterSpectrum.get()[i][1];
            ySpectrum.get()[i][1] = ySpectrum.get()[i][0] * filterSpectrum.get()[i][1] + ySpectrum.get()[i][1] * filterSpectrum.get()[i][0];
            ySpectrum.get()[i][0] = tmp;
        }
    }

    public static void getBestF0Contour(int f0Length, final double[][] f0Candidates, final double[][] f0Scores, int numberOfBands, double[] bestF0Contour) {
        double tmp;
        for (int i = 0; i < f0Length; ++i) {
            tmp = f0Scores[0][i];
            bestF0Contour[i] = f0Candidates[0][i];
            for (int j = 1; j < numberOfBands; ++j) {
                if (tmp > f0Scores[j][i]) {
                    tmp = f0Scores[j][i];
                    bestF0Contour[i] = f0Candidates[j][i];
                }
            }
        }
    }

    public static void fixStep1(final double[] bestF0Contour, int voiceRangeMinimum, double allowedRange, double[] f0Step1) {
        double[] f0Base = new double[bestF0Contour.length];
        for (int i = 0; i < voiceRangeMinimum; ++i)
            f0Base[i] = 0.0;
        if (bestF0Contour.length - voiceRangeMinimum - voiceRangeMinimum >= 0)
            System.arraycopy(bestF0Contour, voiceRangeMinimum, f0Base, voiceRangeMinimum, bestF0Contour.length - voiceRangeMinimum - voiceRangeMinimum);
        for (int i = bestF0Contour.length - voiceRangeMinimum; i < bestF0Contour.length; ++i)
            f0Base[i] = 0.0;

        for (int i = 0; i < voiceRangeMinimum; ++i)
            f0Step1[i] = 0.0;
        for (int i = voiceRangeMinimum; i < bestF0Contour.length; ++i)
            f0Step1[i] = Math.abs((f0Base[i] - f0Base[i - 1]) / (World.MY_SAFE_GUARD_MINIMUM + f0Base[i])) < allowedRange ? f0Base[i] : 0.0;
    }

    public static void fixStep2(final double[] f0Step1, int voiceRangeMinimum, double[] f0Step2) {
        System.arraycopy(f0Step1, 0, f0Step2, 0, f0Step1.length);

        int center = (voiceRangeMinimum - 1) / 2;
        for (int i = center; i < f0Step1.length - center; ++i) {
            for (int j = -center; j <= center; ++j) {
                if (f0Step1[i + j] == 0) {
                    f0Step2[i] = 0.0;
                    break;
                }
            }
        }
    }

    public static void getNumberOfVoicedSections(final double[] f0, int[] positiveIndex, int[] negativeIndex, Holder<Integer> positiveCount, Holder<Integer> negativeCount) {
        positiveCount.set(0);
        negativeCount.set(0);
        for (int i = 1; i < f0.length; ++i)
            if (f0[i] == 0 && f0[i - 1] != 0) {
                negativeIndex[negativeCount.get()] = i - 1;
                negativeCount.set(negativeCount.get() + 1);
            } else {
                if (f0[i - 1] == 0 && f0[i] != 0) {
                    positiveIndex[positiveCount.get()] = i;
                    positiveCount.set(positiveCount.get() + 1);
                }
            }
    }

    public static double selectBestF0(double currentF0, double pastF0, final double[][] f0Candidates, int numberOfCandidates, int targetIndex, double allowedRange) {
        double referenceF0 = (currentF0 * 3.0 - pastF0) / 2.0;

        double minimumError = Math.abs(referenceF0 - f0Candidates[0][targetIndex]);
        double bestF0 = f0Candidates[0][targetIndex];

        double currentError;
        for (int i = 1; i < numberOfCandidates; ++i) {
            currentError = Math.abs(referenceF0 - f0Candidates[i][targetIndex]);
            if (currentError < minimumError) {
                minimumError = currentError;
                bestF0 = f0Candidates[i][targetIndex];
            }
        }
        if (Math.abs(1.0 - bestF0 / referenceF0) > allowedRange)
            return 0.0;
        return bestF0;
    }

    public static void fixStep3(final double[] f0Step2, final double[][] f0Candidates, int numberOfCandidates, double allowedRange, final int[] negativeIndex, int negativeCount, double[] f0Step3) {
        System.arraycopy(f0Step2, 0, f0Step3, 0, f0Step2.length);

        int limit;
        for (int i = 0; i < negativeCount; ++i) {
            limit = i == negativeCount - 1 ? f0Step2.length - 1 : negativeIndex[i + 1];
            for (int j = negativeIndex[i]; j < limit; ++j) {
                f0Step3[j + 1] = selectBestF0(f0Step3[j], f0Step3[j - 1], f0Candidates, numberOfCandidates, j + 1, allowedRange);
                if (f0Step3[j + 1] == 0)
                    break;
            }
        }
    }

    public static void fixStep4(final double[] f0Step3, final double[][] f0Candidates, int numberOfCandidates, double allowedRange, final int[] positiveIndex, int positiveCount, double[] f0Step4) {
        System.arraycopy(f0Step3, 0, f0Step4, 0, f0Step3.length);

        int limit;
        for (int i = positiveCount - 1; i >= 0; --i) {
            limit = i == 0 ? 1 : positiveIndex[i - 1];
            for (int j = positiveIndex[i]; j > limit; --j) {
                f0Step4[j - 1] = selectBestF0(f0Step4[j], f0Step4[j + 1], f0Candidates, numberOfCandidates, j - 1, allowedRange);
                if (f0Step4[j - 1] == 0)
                    break;
            }
        }
    }

    public static void fixF0Contour(double framePeriod, int numberOfCandidates, final double[][] f0Candidates, final double[] bestF0Contour, int f0Length, double f0Floor, double allowedRange, double[] fixedF0Contour) {
        int voiceRangeMinimum = (int) (0.5 + 1000.0 / framePeriod / f0Floor) * 2 + 1;

        if (f0Length <= voiceRangeMinimum) return;

        double[] f0Tmp1 = new double[f0Length];
        double[] f0Tmp2 = new double[f0Length];

        fixStep1(bestF0Contour, voiceRangeMinimum, allowedRange, f0Tmp1);
        fixStep2(f0Tmp1, voiceRangeMinimum, f0Tmp2);

        Holder<Integer> positiveCount = new Holder<>(0), negativeCount = new Holder<>(0);
        int[] positiveIndex = new int[f0Length];
        int[] negativeIndex = new int[f0Length];
        getNumberOfVoicedSections(f0Tmp2, positiveIndex, negativeIndex, positiveCount, negativeCount);
        fixStep3(f0Tmp2, f0Candidates, numberOfCandidates, allowedRange, negativeIndex, negativeCount.get(), f0Tmp1);
        fixStep4(f0Tmp1, f0Candidates, numberOfCandidates, allowedRange, positiveIndex, positiveCount.get(), fixedF0Contour);
    }

    public static void getFilteredSignal(World world, int halfAverageLength, int fftSize, final FFT.FFTComplex ySpectrum, int yLength, double[] filteredSignal) {
        double[] lowPassFilter = new double[fftSize];
        Utils.nuttallWindow(lowPassFilter);
        for (int i = halfAverageLength * 4; i < fftSize; ++i)
            lowPassFilter[i] = 0.0;

        FFT.FFTComplex lowPassFilterSpectrum = new FFT.FFTComplex(fftSize);
        FFT.FFTPlan forwardFFT = FFT.r2c(fftSize, lowPassFilter, lowPassFilterSpectrum, FFT.ESTIMATE);
        world.getFft().execute(forwardFFT);

        double tmp = ySpectrum.get()[0][0] * lowPassFilterSpectrum.get()[0][0] - ySpectrum.get()[0][1] * lowPassFilterSpectrum.get()[0][1];
        lowPassFilterSpectrum.get()[0][1] = ySpectrum.get()[0][0] * lowPassFilterSpectrum.get()[0][1] + ySpectrum.get()[0][1] * lowPassFilterSpectrum.get()[0][0];
        lowPassFilterSpectrum.get()[0][0] = tmp;
        for (int i = 1; i <= fftSize / 2; ++i) {
            tmp = ySpectrum.get()[i][0] * lowPassFilterSpectrum.get()[i][0] - ySpectrum.get()[i][1] * lowPassFilterSpectrum.get()[i][1];
            lowPassFilterSpectrum.get()[i][1] = ySpectrum.get()[i][0] * lowPassFilterSpectrum.get()[i][1] + ySpectrum.get()[i][1] * lowPassFilterSpectrum.get()[i][0];
            lowPassFilterSpectrum.get()[i][0] = tmp;
            lowPassFilterSpectrum.get()[fftSize - i - 1][0] = lowPassFilterSpectrum.get()[i][0];
            lowPassFilterSpectrum.get()[fftSize - i - 1][1] = lowPassFilterSpectrum.get()[i][1];
        }

        FFT.FFTPlan inverseFFT = FFT.c2r(fftSize, lowPassFilterSpectrum, filteredSignal, FFT.ESTIMATE);
        world.getFft().execute(inverseFFT);

        int indexBias = halfAverageLength * 2;
        if (yLength >= 0) System.arraycopy(filteredSignal, indexBias, filteredSignal, 0, yLength);
    }

    public static int checkEvent(int x) {
        return x > 0 ? 1 : 0;
    }

    public static int zeroCrossingEngine(final double[] filteredSignal, int yLength, double fs, double[] intervalLocations, double[] intervals) {
        int[] negativeGoingPoints = new int[yLength];

        for (int i = 0; i < yLength - 1; ++i)
            negativeGoingPoints[i] = 0.0 < filteredSignal[i] && filteredSignal[i + 1] <= 0.0 ? i + 1 : 0;
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
            fineEdges[i] = edges[i] - filteredSignal[edges[i] - 1] / (filteredSignal[edges[i]] - filteredSignal[edges[i] - 1]);

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

        for (int i = 0; i < yLength - 1; ++i)
            filteredSignal[i] = filteredSignal[i] - filteredSignal[i + 1];
        zeroCrossings.setNumberOfPeaks(zeroCrossingEngine(filteredSignal, yLength - 1, actualFs, zeroCrossings.getPeakIntervalLocations(), zeroCrossings.getPeakIntervals()));

        for (int i = 0; i < yLength - 1; ++i)
            filteredSignal[i] = -filteredSignal[i];
        zeroCrossings.setNumberOfDips(zeroCrossingEngine(filteredSignal, yLength - 1, actualFs, zeroCrossings.getDipIntervalLocations(), zeroCrossings.getDipIntervals()));
    }

    public static void getF0CandidateContourSub(final double[][] interpolatedF0Set, int f0Length, double f0Floor, double f0Ceil, double boundaryF0, double[] f0Candidate, double[] f0Score) {
        for (int i = 0; i < f0Length; ++i) {
            f0Candidate[i] = (interpolatedF0Set[0][i] + interpolatedF0Set[1][i] + interpolatedF0Set[2][i] + interpolatedF0Set[3][i]) / 4.0;

            f0Score[i] = Math.sqrt(
                    ((interpolatedF0Set[0][i] - f0Candidate[i]) * (interpolatedF0Set[0][i] - f0Candidate[i])
                            + (interpolatedF0Set[1][i] - f0Candidate[i]) * (interpolatedF0Set[1][i] - f0Candidate[i])
                            + (interpolatedF0Set[2][i] - f0Candidate[i]) * (interpolatedF0Set[2][i] - f0Candidate[i])
                            + (interpolatedF0Set[3][i] - f0Candidate[i]) * (interpolatedF0Set[3][i] - f0Candidate[i])) / 3.0);

            if (f0Candidate[i] > boundaryF0 || f0Candidate[i] < boundaryF0 / 2.0 || f0Candidate[i] > f0Ceil || f0Candidate[i] < f0Floor) {
                f0Candidate[i] = 0.0;
                f0Score[i] = World.MAXIMUM_VALUE;
            }
        }
    }

    public static void getF0CandidateContour(final ZeroCrossings zeroCrossings, double boundaryF0, double f0Floor, double f0Ceil, final double[] temporalPositions, int f0Length, double[] f0Candidate, double[] f0Score) {
        if (0 == checkEvent(zeroCrossings.getNumberOfNegatives() - 2) * checkEvent(zeroCrossings.getNumberOfPositives() - 2)
                * checkEvent(zeroCrossings.getNumberOfPeaks() - 2) * checkEvent(zeroCrossings.getNumberOfDips() - 2)) {
            for (int i = 0; i < f0Length; ++i) {
                f0Score[i] = World.MAXIMUM_VALUE;
                f0Candidate[i] = 0.0;
            }
            return;
        }

        double[][] interpolatedF0Set = new double[f0Length][4];
        for (int i = 0; i < 4; ++i)
            interpolatedF0Set[i] = new double[f0Length];

        Utils.interp1(zeroCrossings.getNegativeIntervalLocations(), zeroCrossings.getNegativeIntervals(), temporalPositions, interpolatedF0Set[0]);
        Utils.interp1(zeroCrossings.getPositiveIntervalLocations(), zeroCrossings.getPositiveIntervals(), temporalPositions, interpolatedF0Set[1]);
        Utils.interp1(zeroCrossings.getPeakIntervalLocations(), zeroCrossings.getPeakIntervals(), temporalPositions, interpolatedF0Set[2]);
        Utils.interp1(zeroCrossings.getDipIntervalLocations(), zeroCrossings.getDipIntervals(), temporalPositions, interpolatedF0Set[3]);

        getF0CandidateContourSub(interpolatedF0Set, f0Length, f0Floor, f0Ceil, boundaryF0, f0Candidate, f0Score);
    }

    public static void getF0CandidateFromRawEvent(World world, double boundaryF0, double fs, final FFT.FFTComplex ySpectrum, int yLength, int fftSize, double f0Floor,
                                                  double f0Ceil, final double[] temporalPositions, int f0Length, double[] f0Score, double[] f0Candidate) {
        double[] filteredSignal = new double[fftSize];
        getFilteredSignal(world, (int) Math.round(fs / boundaryF0 / 2.0), fftSize, ySpectrum, yLength, filteredSignal);

        ZeroCrossings zeroCrossings = new ZeroCrossings();
        getFourZeroCrossingIntervals(filteredSignal, yLength, fs, zeroCrossings);

        getF0CandidateContour(zeroCrossings, boundaryF0, f0Floor, f0Ceil, temporalPositions, f0Length, f0Candidate, f0Score);
    }

    public static void getF0CandidatesAndScores(World world, final double[] boundaryF0List, int numberOfBands, double actualFs, int yLength, final double[] temporalPositions, int f0Length,
                                                final FFT.FFTComplex ySpectrum, int fftSize, double f0Floor, double f0Ceil, double[][] rawF0Candidates, double[][] rawF0Scores) {
        double[] f0Candidate = new double[f0Length];
        double[] f0Score = new double[f0Length];

        for (int i = 0; i < numberOfBands; ++i) {
            getF0CandidateFromRawEvent(world, boundaryF0List[i], actualFs, ySpectrum, yLength, fftSize, f0Floor, f0Ceil, temporalPositions, f0Length, f0Score, f0Candidate);
            for (int j = 0; j < f0Length; ++j) {
                rawF0Scores[i][j] = f0Score[j] / (f0Candidate[j] + World.MY_SAFE_GUARD_MINIMUM);
                rawF0Candidates[i][j] = f0Candidate[j];
            }
        }

    }

    public static void dioGeneralBody(World world, final double[] x, int fs, double framePeriod, double f0Floor, double f0Ceil, double channelsInOctave, int speed, double allowedRange, double[] temporalPositions, double[] f0) {
        int numberOfBands = 1 + (int) (Math.log(f0Ceil / f0Floor) / World.LOG_2 * channelsInOctave);
        double[] boundaryF0List = new double[numberOfBands];
        for (int i = 0; i < numberOfBands; ++i)
            boundaryF0List[i] = f0Floor * Math.pow(2.0, (i + 1) / channelsInOctave);

        int decimationRatio = Math.max(Math.min(speed, 12), 1);
        int yLength = (1 + (x.length / decimationRatio));
        double actualFs = (double) (fs) / decimationRatio;
        int fftSize = Utils.getSuitableFFTSize((int) (yLength + Math.round(actualFs / World.CUT_OFF) * 2 + 1 + (4 * (int) (1.0 + actualFs / boundaryF0List[0] / 2.0))));

        FFT.FFTComplex ySpectrum = new FFT.FFTComplex(fftSize);
        getSpectrumForEstimation(world, x, yLength, actualFs, fftSize, decimationRatio, ySpectrum);

        int f0Length = getSamplesForDIO(fs, x.length, framePeriod);
        double[][] f0Candidates = new double[f0Length][numberOfBands];
        double[][] f0Scores = new double[f0Length][numberOfBands];
        for (int i = 0; i < numberOfBands; ++i) {
            f0Candidates[i] = new double[f0Length];
            f0Scores[i] = new double[f0Length];
        }

        for (int i = 0; i < f0Length; ++i)
            temporalPositions[i] = i * framePeriod / 1000.0;

        getF0CandidatesAndScores(world, boundaryF0List, numberOfBands, actualFs, yLength, temporalPositions, f0Length, ySpectrum, fftSize, f0Floor, f0Ceil, f0Candidates, f0Scores);

        double[] bestF0Contour = new double[f0Length];
        getBestF0Contour(f0Length, f0Candidates, f0Scores, numberOfBands, bestF0Contour);

        fixF0Contour(framePeriod, numberOfBands, f0Candidates, bestF0Contour, f0Length, f0Floor, allowedRange, f0);
    }

    public static int getSamplesForDIO(int fs, int xLength, double framePeriod) {
        return (int) (1000.0 * xLength / fs / framePeriod) + 1;
    }

    public void dio(final double[] x, int fs, final DioOption option, double[] temporalPositions, double[] f0) {
        dioGeneralBody(world, x, fs, option.framePeriod, option.f0Floor, option.f0Ceil, option.channelsInOctave, option.speed, option.allowedRange, temporalPositions, f0);
    }

    public void initializeDioOption(DioOption option) {
        option.channelsInOctave = 2.0;
        option.f0Ceil = World.CEIL_F0;
        option.f0Floor = World.FLOOR_F0;
        option.framePeriod = 5;

        option.speed = 1;
        option.allowedRange = 0.1;
    }

    @Data
    @NoArgsConstructor
    @AllArgsConstructor
    public static class DioOption {
        double f0Floor;
        double f0Ceil;
        double channelsInOctave;
        double framePeriod;
        int speed;
        double allowedRange;
    }
}
