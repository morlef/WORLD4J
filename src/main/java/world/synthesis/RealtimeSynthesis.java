package world.synthesis;

import world.World;
import world.common.ForwardRealFFT;
import world.common.InverseRealFFT;
import world.common.MinimumPhaseAnalysis;
import world.util.Holder;
import world.util.Utils;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/synthrealtime.cpp
 */
public record RealtimeSynthesis(World world) {
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

    public static void getAperiodicResponse(World world, int noiseSize, int fftSize, final double[] spectrum, final double[] aperiodicRatio, double currentVuv,
                                     final ForwardRealFFT forwardRealFFT, final InverseRealFFT inverseRealFFT, final MinimumPhaseAnalysis minimumPhase, double[] aperiodicResponse) {
        getNoiseSpectrum(world, noiseSize, fftSize, forwardRealFFT);

        if (currentVuv != 0.0)
            for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i)
                minimumPhase.getLogSpectrum()[i] = Math.log(spectrum[i] * aperiodicRatio[i] + World.MY_SAFE_GUARD_MINIMUM) / 2.0;
        else
            for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i)
                minimumPhase.getLogSpectrum()[i] = Math.log(spectrum[i]) / 2.0;
        Utils.getMinimumPhaseSpectrum(world, minimumPhase);

        for (int i = 0; i <= fftSize / 2; ++i) {
            inverseRealFFT.getSpectrum().get()[i][0] =
                    minimumPhase.getMinimumPhaseSpectrum().get()[i][0] *
                            forwardRealFFT.getSpectrum().get()[i][0] -
                                    minimumPhase.getMinimumPhaseSpectrum().get()[i][1] *
                                            forwardRealFFT.getSpectrum().get()[i][1];
            inverseRealFFT.getSpectrum().get()[i][1] =
                    minimumPhase.getMinimumPhaseSpectrum().get()[i][0] *
                            forwardRealFFT.getSpectrum().get()[i][1] +
                                    minimumPhase.getMinimumPhaseSpectrum().get()[i][1] *
                                            forwardRealFFT.getSpectrum().get()[i][0];
        }
        world.getFft().execute(inverseRealFFT.getInverseFFT());
        Utils.fftShift(inverseRealFFT.getWaveform(), aperiodicResponse);
    }

    public static void clearRingBuffer(int start, int end, WorldSynthesizer synth) {
        int pointer;
        for (int i = start; i < end; ++i) {
            pointer = i % synth.numberOfPointers;
            synth.numberOfPulses[pointer] = 0;
            if (synth.pulseLocations[pointer] != null) {
                synth.pulseLocations[pointer] = null;
            }
            if (synth.interpolatedVuv[pointer] != null) {
                synth.interpolatedVuv[pointer] = null;
            }
            if (synth.pulseLocationsIndex[pointer] != null) {
                synth.pulseLocationsIndex[pointer] = null;
            }
        }
    }

    public static int seekSynthesizer(double currentLocation, WorldSynthesizer synth) {
        int frameNumber = (int) (currentLocation / synth.framePeriod);

        int tmpPointer = synth.currentPointer2;
        int tmp;
        for (int i = 0; i < synth.headPointer - synth.currentPointer2; ++i) {
            tmp = (tmpPointer + i) % synth.numberOfPointers;
            if (synth.f0Origin[tmp] <= frameNumber &&
                    frameNumber < synth.f0Origin[tmp] + synth.f0Length[tmp]) {
                tmpPointer += i;
                break;
            }
        }
        clearRingBuffer(synth.currentPointer2, tmpPointer, synth);
        synth.currentPointer2 = tmpPointer;
        return 1;
    }

    public static void searchPointer(int frame, WorldSynthesizer synth, int flag, Holder<double[]> front, Holder<double[]> next) {
        int pointer = synth.currentPointer2 % synth.numberOfPointers;
        int index = -1;
        for (int i = 0; i < synth.f0Length[pointer]; ++i)
            if (synth.f0Origin[pointer] + i == frame) {
                index = i;
                break;
            }

        double[][][] tmpPointer = flag == 0 ? synth.spectrogram : synth.aperiodicity;

        front.set(tmpPointer[pointer][index]);
        next.set(index == synth.f0Length[pointer] - 1 ? tmpPointer[(synth.currentPointer2 + 1) % synth.numberOfPointers][0] : tmpPointer[pointer][index + 1]);
    }

    public static void getPeriodicResponse(World world, int fftSize, final double[] spectrum, final double[] aperiodicRatio, double currentVuv, final InverseRealFFT inverseRealFFT, final MinimumPhaseAnalysis minimumPhase, final double[] periodicResponse) {
        if (currentVuv <= 0.5 || aperiodicRatio[0] > 0.999) {
            for (int i = 0; i < fftSize; ++i) periodicResponse[i] = 0.0;
            return;
        }

        for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i)
            minimumPhase.getLogSpectrum()[i] = Math.log(spectrum[i] * (1.0 - aperiodicRatio[i]) + World.MY_SAFE_GUARD_MINIMUM) / 2.0;
        Utils.getMinimumPhaseSpectrum(world, minimumPhase);

        for (int i = 0; i <= fftSize / 2; ++i) {
            inverseRealFFT.getSpectrum().get()[i][0] =
                    minimumPhase.getMinimumPhaseSpectrum().get()[i][0];
            inverseRealFFT.getSpectrum().get()[i][1] =
                    minimumPhase.getMinimumPhaseSpectrum().get()[i][1];
        }

        world.getFft().execute(inverseRealFFT.getInverseFFT());
        Utils.fftShift(inverseRealFFT.getWaveform(), periodicResponse);
    }

    public static void getSpectralEnvelope(double currentLocation, WorldSynthesizer synth, double[] spectralEnvelope) {
        int currentFrameFloor = (int) (currentLocation / synth.framePeriod);

        int currentFrameCeil = (int) (Math.ceil(currentLocation / synth.framePeriod));
        double interpolation = currentLocation / synth.framePeriod - currentFrameFloor;
        Holder<double[]> front = new Holder<>();
        Holder<double[]> next = new Holder<>();
        searchPointer(currentFrameFloor, synth, 0, front, next);

        if (currentFrameFloor == currentFrameCeil)
            for (int i = 0; i <= synth.fftSize / 2; ++i)
                spectralEnvelope[i] = Math.abs(front.get()[i]);
        else
            for (int i = 0; i <= synth.fftSize / 2; ++i)
                spectralEnvelope[i] =
                        (1.0 - interpolation) * Math.abs(front.get()[i]) + interpolation * Math.abs(next.get()[i]);
    }

    public static void getAperiodicRatio(double currentLocation, WorldSynthesizer synth, double[] aperiodicSpectrum) {
        int currentFrameFloor = (int) (currentLocation / synth.framePeriod);

        int currentFrameCeil = (int) (Math.ceil(currentLocation / synth.framePeriod));
        double interpolation = currentLocation / synth.framePeriod - currentFrameFloor;

        Holder<double[]> front = new Holder<>();
        Holder<double[]> next = new Holder<>();
        searchPointer(currentFrameFloor, synth, 1, front, next);

        if (currentFrameFloor == currentFrameCeil)
            for (int i = 0; i <= synth.fftSize / 2; ++i)
                aperiodicSpectrum[i] = Math.pow(Utils.getSafeAperiodicity(front.get()[i]), 2.0);
        else
            for (int i = 0; i <= synth.fftSize / 2; ++i)
                aperiodicSpectrum[i] = Math.pow((1.0 - interpolation) * Utils.getSafeAperiodicity(front.get()[i]) + interpolation * Utils.getSafeAperiodicity(next.get()[i]), 2.0);
    }

    public static double getCurrentVUV(int currentLocation, WorldSynthesizer synth) {
        double currentVuv;
        int pointer = synth.currentPointer % synth.numberOfPointers;

        int startSample = Math.max(0, (int) (Math.ceil((synth.f0Origin[pointer] - 1) * synth.framePeriod * synth.fs)));

        currentVuv = synth.interpolatedVuv[pointer][currentLocation - startSample + 1];
        return currentVuv;
    }

    public static void getOneFrameSegment(World world, int noiseSize, int currentLocation, WorldSynthesizer synth) {
        double[] aperiodicResponse = new double[synth.fftSize];
        double[] periodicResponse = new double[synth.fftSize];
        double[] spectralEnvelope = new double[synth.fftSize];
        double[] aperiodicRatio = new double[synth.fftSize];

        double tmpLocation = (double) (currentLocation) / synth.fs;
        seekSynthesizer(tmpLocation, synth);
        getSpectralEnvelope(tmpLocation, synth, spectralEnvelope);
        getAperiodicRatio(tmpLocation, synth, aperiodicRatio);

        double currentVuv = getCurrentVUV(currentLocation, synth);

        getPeriodicResponse(world, synth.fftSize, spectralEnvelope, aperiodicRatio, currentVuv, synth.inverseRealFFT, synth.minimumPhase, periodicResponse);

        getAperiodicResponse(world, noiseSize, synth.fftSize, spectralEnvelope, aperiodicRatio, currentVuv, synth.forwardRealFFT, synth.inverseRealFFT, synth.minimumPhase, aperiodicResponse);

        double sqrtNoiseSize = Math.sqrt(noiseSize);
        for (int i = 0; i < synth.fftSize; ++i)
            synth.impulseResponse[i] = (periodicResponse[i] * sqrtNoiseSize + aperiodicResponse[i]) / synth.fftSize;
    }

    public static void getTemporalParametersForTimeBase(final double[] f0, int f0Length, WorldSynthesizer synth, double[] coarseTimeAxis, double[] coarseF0, double[] coarseVuv) {
        int cumulativeFrame = Math.max(0, synth.cumulativeFrame - f0Length);
        coarseF0[0] = synth.handoffF0;
        coarseTimeAxis[0] = cumulativeFrame * synth.framePeriod;
        coarseVuv[0] = synth.handoffF0 == 0 ? 0.0 : 1.0;
        for (int i = 0; i < f0Length; ++i) {
            coarseTimeAxis[i + synth.handoff] = (i + cumulativeFrame + synth.handoff) * synth.framePeriod;
            coarseF0[i + synth.handoff] = f0[i];
            coarseVuv[i + synth.handoff] = f0[i] == 0.0 ? 0.0 : 1.0;
        }
    }

    public static void getPulseLocationsForTimeBase(final double[] interpolatedF0, final double[] timeAxis, int numberOfSamples, WorldSynthesizer synth) {
        double[] totalPhase = new double[numberOfSamples + synth.handoff];
        totalPhase[0] = synth.handoff == 1 ? synth.handoffPhase : 2.0 * Math.PI * interpolatedF0[0] / synth.fs;

        totalPhase[1] = totalPhase[0] + 2.0 * Math.PI * interpolatedF0[0] / synth.fs;
        for (int i = 1 + synth.handoff; i < numberOfSamples + synth.handoff; ++i)
            totalPhase[i] = totalPhase[i - 1] + 2.0 * Math.PI * interpolatedF0[i - synth.handoff] / synth.fs;
        synth.handoffPhase = totalPhase[numberOfSamples - 1 + synth.handoff];

        double[] wrapPhase = new double[numberOfSamples + synth.handoff];
        for (int i = 0; i < numberOfSamples + synth.handoff; ++i)
            wrapPhase[i] = (totalPhase[i] % 2.0 * Math.PI);

        double[] wrapPhaseAbs = new double[numberOfSamples + synth.handoff];
        for (int i = 0; i < numberOfSamples - 1 + synth.handoff; ++i)
            wrapPhaseAbs[i] = Math.abs(wrapPhase[i + 1] - wrapPhase[i]);

        int pointer = synth.headPointer % synth.numberOfPointers;
        int numberOfPulses = 0;
        for (int i = 0; i < numberOfSamples - 1 + synth.handoff; ++i)
            if (wrapPhaseAbs[i] > Math.PI) {
                synth.pulseLocations[pointer][numberOfPulses] = timeAxis[i] - (double) (synth.handoff) / synth.fs;
                synth.pulseLocationsIndex[pointer][numberOfPulses] = (int) Math.round(synth.pulseLocations[pointer][numberOfPulses] * synth.fs);
                ++numberOfPulses;
            }
        synth.numberOfPulses[pointer] = numberOfPulses;

        if (numberOfPulses != 0)
            synth.lastLocation = synth.pulseLocationsIndex[pointer][numberOfPulses - 1];
    }

    public static void getTimeBase(final double[] f0, int f0Length, int startSample, int numberOfSamples, WorldSynthesizer synth) {
        double[] coarseTimeAxis = new double[f0Length + synth.handoff];
        double[] coarseF0 = new double[f0Length + synth.handoff];
        double[] coarseVuv = new double[f0Length + synth.handoff];

        getTemporalParametersForTimeBase(f0, f0Length, synth, coarseTimeAxis, coarseF0, coarseVuv);

        double[] interpolatedF0 = new double[numberOfSamples];
        double[] timeAxis = new double[numberOfSamples];

        for (int i = 0; i < numberOfSamples; ++i)
            timeAxis[i] = (i + startSample) / (double) (synth.fs);

        int pointer = synth.headPointer % synth.numberOfPointers;
        Utils.interp1(coarseTimeAxis, coarseF0, timeAxis, interpolatedF0);
        Utils.interp1(coarseTimeAxis, coarseVuv, timeAxis, synth.interpolatedVuv[pointer]);
        for (int i = 0; i < numberOfSamples; ++i) {
            synth.interpolatedVuv[pointer][i] = synth.interpolatedVuv[pointer][i] > 0.5 ? 1.0 : 0.0;
            interpolatedF0[i] = synth.interpolatedVuv[pointer][i] == 0.0 ? World.DEFAULT_F0 : interpolatedF0[i];
        }

        getPulseLocationsForTimeBase(interpolatedF0, timeAxis, numberOfSamples, synth);

        synth.handoffF0 = interpolatedF0[numberOfSamples - 1];
    }

    public static int getNextPulseLocationIndex(WorldSynthesizer synth) {
        int pointer = synth.currentPointer % synth.numberOfPointers;
        if (synth.i < synth.numberOfPulses[pointer] - 1)
            return synth.pulseLocationsIndex[pointer][synth.i + 1];
        else if (synth.currentPointer == synth.headPointer - 1)
            return 0;

        for (int i = 1; i < synth.numberOfPointers; ++i) {
            pointer = (i + synth.currentPointer) % synth.numberOfPointers;
            if (synth.numberOfPulses[pointer] != 0)
                return synth.pulseLocationsIndex[pointer][0];
        }
        return 0;
    }

    public static int updateSynthesizer(WorldSynthesizer synth) {
        int pointer = synth.currentPointer % synth.numberOfPointers;
        if (synth.i < synth.numberOfPulses[pointer] - 1) {
            synth.i++;
            return 1;
        } else {
            if (synth.currentPointer == synth.headPointer - 1) return 0;
        }

        for (int i = 1; i < synth.numberOfPointers; ++i) {
            pointer = (i + synth.currentPointer) % synth.numberOfPointers;
            if (synth.numberOfPulses[pointer] != 0) {
                synth.i = 0;
                synth.currentPointer += i;
                return 1;
            }
        }
        return 0;
    }

    public static int checkSynthesizer(WorldSynthesizer synth) {
        if (synth.synthesizedSample + synth.bufferSize >= synth.lastLocation)
            return 0;

        int pointer = synth.currentPointer % synth.numberOfPointers;
        while (synth.numberOfPulses[pointer] == 0) {
            if (synth.currentPointer == synth.headPointer) break;
            synth.currentPointer++;
            pointer = synth.currentPointer % synth.numberOfPointers;
        }
        return 1;
    }

    public int addParameters(double[] f0, int f0Length, double[][] spectrogram, double[][] aperiodicity, WorldSynthesizer synth) {
        if (synth.headPointer - synth.currentPointer2 == synth.numberOfPointers)
            return 0;
        int pointer = synth.headPointer % synth.numberOfPointers;
        synth.f0Length[pointer] = f0Length;
        synth.f0Origin[pointer] = synth.cumulativeFrame + 1;
        synth.cumulativeFrame += f0Length;

        synth.spectrogram[pointer] = spectrogram;
        synth.aperiodicity[pointer] = aperiodicity;
        if (synth.cumulativeFrame < 1) {
            synth.handoffF0 = f0[f0Length - 1];
            synth.numberOfPulses[pointer] = 0;
            synth.headPointer++;
            synth.handoff = 1;
            return 1;
        }

        int startSample = Math.max(0, (int) (Math.ceil((synth.cumulativeFrame - f0Length) * synth.framePeriod * synth.fs)));
        int endSample = (int) (Math.ceil((synth.cumulativeFrame) * synth.framePeriod * synth.fs));
        int numberOfSamples = endSample - startSample;

        synth.interpolatedVuv[pointer] = new double[numberOfSamples + 1];
        synth.pulseLocations[pointer] = new double[numberOfSamples];
        synth.pulseLocationsIndex[pointer] = new int[numberOfSamples];

        getTimeBase(f0, f0Length, startSample, numberOfSamples, synth);

        synth.handoffF0 = f0[f0Length - 1];
        synth.headPointer++;
        synth.handoff = 1;
        return 1;
    }

    public int isLocked(WorldSynthesizer synth) {
        int judge = 0;
        if (synth.headPointer - synth.currentPointer2 == synth.numberOfPointers)
            judge++;
        if (synth.synthesizedSample + synth.bufferSize >= synth.lastLocation)
            judge++;

        return judge == 2 ? 1 : 0;
    }

    public int synthesis2(WorldSynthesizer synth) {
        if (checkSynthesizer(synth) == 0)
            return 0;
        if (synth.bufferSize + synth.fftSize >= 0)
            System.arraycopy(synth.buffer, synth.bufferSize, synth.buffer, 0, synth.bufferSize + synth.fftSize);

        int pointer = synth.currentPointer % synth.numberOfPointers;
        int noiseSize, offset, tmp, index;
        int currentLocation = synth.pulseLocationsIndex[pointer][synth.i];
        while (currentLocation < synth.synthesizedSample + synth.bufferSize) {
            tmp = getNextPulseLocationIndex(synth);
            noiseSize = tmp - currentLocation;

            getOneFrameSegment(world, noiseSize, currentLocation, synth);
            offset =
                    currentLocation - synth.synthesizedSample - synth.fftSize / 2 + 1;
            for (int i = Math.max(0, -offset); i < synth.fftSize; ++i) {
                index = i + offset;
                synth.buffer[index] += synth.impulseResponse[i];
            }
            currentLocation = tmp;
            updateSynthesizer(synth);
        }
        synth.synthesizedSample += synth.bufferSize;
        seekSynthesizer(synth.synthesizedSample, synth);
        return 1;
    }
}
