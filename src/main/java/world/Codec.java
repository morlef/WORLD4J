package world;

import lombok.experimental.PackagePrivate;
import world.common.ForwardRealFFT;
import world.common.InverseComplexFFT;
import world.util.Utils;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/codec.cpp
 */
@PackagePrivate
public record Codec(World world) {
    public int getNumberOfAperiodicities(int fs) {
        return  (int)(Math.min(World.UPPER_LIMIT, fs / 2.0 - World.FREQUENCY_INTERVAL) / World.FREQUENCY_INTERVAL);
    }

    public void codeAperiodicity(final double[][] aperiodicity, int f0Length, int fs, int fftSize, double[][] codedAperiodicity) {
        int numberOfAperiodicities = getNumberOfAperiodicities(fs);
        double[] coarseFrequencyAxis = new double[numberOfAperiodicities];
        for (int i = 0; i < numberOfAperiodicities; ++i)
            coarseFrequencyAxis[i] = World.FREQUENCY_INTERVAL * (i + 1.0);

        double[] logAperiodicity = new double[fftSize / 2 + 1];

        for (int i = 0; i < f0Length; ++i) {
            for (int j = 0; j < fftSize / 2 + 1; ++j)
                logAperiodicity[j] = 20 * Math.log10(aperiodicity[i][j]);
            Utils.interp1Q(0, (double) (fs) / fftSize, logAperiodicity, coarseFrequencyAxis, codedAperiodicity[i]);
        }
    }

    public void decodeAperiodicity(final double[][] codedAperiodicity, int f0Length, int fs, int fftSize, double[][] aperiodicity) {
        initializeAperiodicity(f0Length, fftSize, aperiodicity);
        int numberOfAperiodicities = getNumberOfAperiodicities(fs);
        double[] frequencyAxis = new double[fftSize / 2 + 1];
        for (int i = 0; i <= fftSize / 2; ++i)
            frequencyAxis[i] = (double) (fs) / fftSize * i;

        double[] coarseFrequencyAxis = new double[numberOfAperiodicities + 2];
        for (int i = 0; i <= numberOfAperiodicities; ++i)
            coarseFrequencyAxis[i] = i * World.FREQUENCY_INTERVAL;
        coarseFrequencyAxis[numberOfAperiodicities + 1] = fs / 2.0;

        double[] coarseAperiodicity = new double[numberOfAperiodicities + 2];
        coarseAperiodicity[0] = -60.0;
        coarseAperiodicity[numberOfAperiodicities + 1] = -World.MY_SAFE_GUARD_MINIMUM;

        for (int i = 0; i < f0Length; ++i) {
            if (checkVUV(codedAperiodicity[i], numberOfAperiodicities, coarseAperiodicity) == 1) continue;
            getAperiodicity(coarseFrequencyAxis, coarseAperiodicity, frequencyAxis, fftSize, aperiodicity[i]);
        }
    }

    public void codeSpectralEnvelope(final double[][] spectrogram, int f0Length, int fs, int fftSize, int numberOfDimensions, double[][] codedSpectralEnvelope) {
        double[] melAxis = new double[fftSize / 2];
        double[] frequencyAxis = new double[fftSize / 2 + 1];
        double[] tmpSpectrum = new double[fftSize / 2 + 1];
        FFT.FFTComplex weight = new FFT.FFTComplex(fftSize / 2);

        getParametersForCoding(World.FLOOR_FREQUENCY, Math.min(fs / 2.0, World.CEIL_FREQUENCY), fs, fftSize, melAxis, frequencyAxis, weight);

        ForwardRealFFT forwardRealFFT = new ForwardRealFFT(fftSize / 2);

        for (int i = 0; i < f0Length; ++i) {
            for (int j = 0; j < fftSize / 2 + 1; ++j)
                tmpSpectrum[j] = Math.log(spectrogram[i][j]);
            codeOneFrame(world, tmpSpectrum, frequencyAxis, melAxis, weight, fftSize / 2, numberOfDimensions, forwardRealFFT, codedSpectralEnvelope[i]);
        }
    }

    public void decodeSpectralEnvelope(final double[][] codedSpectralEnvelope, int f0Length, int fs, int fftSize, int numberOfDimensions, double[][] spectrogram) {
        double[] melAxis = new double[fftSize / 2 + 2];
        double[] frequencyAxis = new double[fftSize / 2 + 1];
        FFT.FFTComplex weight = new FFT.FFTComplex(fftSize / 2);

        getParametersForDecoding(World.FLOOR_FREQUENCY, Math.min(fs / 2.0, World.CEIL_FREQUENCY), fs, fftSize, numberOfDimensions, melAxis, frequencyAxis, weight);

        InverseComplexFFT inverseComplexFFT = new InverseComplexFFT(fftSize / 2);

        for (int i = 0; i < f0Length; ++i) {
            decodeOneFrame(world, codedSpectralEnvelope[i], frequencyAxis, fftSize, melAxis, weight, fftSize / 2, numberOfDimensions, inverseComplexFFT, spectrogram[i]);
        }
    }

    public static void initializeAperiodicity(int f0Length, int fftSize, double[][] aperiodicity) {
        for (int i = 0; i < f0Length; ++i)
            for (int j = 0; j < fftSize / 2 + 1; ++j)
                aperiodicity[i][j] = 1.0 - World.MY_SAFE_GUARD_MINIMUM;
    }

    public static int checkVUV(final double[] coarseAperiodicity, int numberOfAperiodicities, double[] tmpAperiodicity) {
        double tmp = 0.0;
        for (int i = 0; i < numberOfAperiodicities; ++i) {
            tmp += coarseAperiodicity[i];
            tmpAperiodicity[i + 1] = coarseAperiodicity[i];
        }
        tmp /= numberOfAperiodicities;

        return tmp > -0.5 ? 1 : 0;
    }

    public static void getAperiodicity(final double[] coarseFrequencyAxis, final double[] coarseAperiodicity, final double[] frequencyAxis, int fftSize, double[] aperiodicity) {
        Utils.interp1(coarseFrequencyAxis, coarseAperiodicity, frequencyAxis, aperiodicity);
        for (int i = 0; i <= fftSize / 2; ++i)
            aperiodicity[i] = Math.pow(10.0, aperiodicity[i] / 20.0);
    }

    public static double frequencyToMel(double frequency) {
        return World.M0 * Math.log(frequency / World.F0 + 1.0);
    }

    public static double melToFrequency(double mel) {
        return World.F0 * (Math.exp(mel / World.M0) - 1.0);
    }

    public static void dctForCodec(World world, final double[] melSpectrum, int maxDimension, final FFT.FFTComplex weight, final ForwardRealFFT forwardRealFFT, int numberOfDimensions, double[] melCepstrum) {
        int bias = maxDimension / 2;
        for (int i = 0; i < maxDimension / 2; ++i) {
            forwardRealFFT.getWaveform()[i] = melSpectrum[i * 2];
            forwardRealFFT.getWaveform()[i + bias] = melSpectrum[maxDimension - (i * 2) - 1];
        }
        world.getFft().execute(forwardRealFFT.getForwardFFT());

        double normalization = Math.sqrt(forwardRealFFT.getFftSize());
        for (int i = 0; i < numberOfDimensions; ++i)
            melCepstrum[i] = (forwardRealFFT.getSpectrum().get()[i][0] * weight.get()[i][0] - forwardRealFFT.getSpectrum().get()[i][1] * weight.get()[i][1]) / normalization;
    }

    public static void iDCTForCodec(World world, final double[] melCepstrum, int maxDimension, final FFT.FFTComplex weight, final InverseComplexFFT inverseComplexFFT, int numberOfDimensions, double[] melSpectrum, int s) {
        double normalization = Math.sqrt(inverseComplexFFT.getFftSize());
        for (int i = 0; i < numberOfDimensions; ++i) {
            inverseComplexFFT.getInput().get()[i][0] = melCepstrum[i + s] * weight.get()[i][0] * normalization;
            inverseComplexFFT.getInput().get()[i][1] = -melCepstrum[i + s] * weight.get()[i][1] * normalization;
        }
        for (int i = numberOfDimensions; i < maxDimension; ++i) {
            inverseComplexFFT.getInput().get()[i][0] = 0.0;
            inverseComplexFFT.getInput().get()[i][1] = 0.0;
        }

        world.getFft().execute(inverseComplexFFT.getInverseFFT());

        for (int i = 0; i < maxDimension / 2; ++i) {
            melSpectrum[i * 2 + s] = inverseComplexFFT.getOutput().get()[i][0];
            melSpectrum[(i * 2) + 1 + s] = inverseComplexFFT.getOutput().get()[maxDimension - i - 1][0];
        }
    }

    public static void codeOneFrame(World world, final double[] logSpectralEnvelope, final double[] frequencyAxis, final double[] melAxis, final FFT.FFTComplex weight, int maxDimension, int numberOfDimensions, final ForwardRealFFT forwardRealFFT, double[] codedSpectralEnvelope) {
        double[] melSpectrum = new double[maxDimension];
        Utils.interp1(frequencyAxis, logSpectralEnvelope, melAxis, melSpectrum);

        dctForCodec(world, melSpectrum, maxDimension, weight, forwardRealFFT, numberOfDimensions, codedSpectralEnvelope);
    }

    public static void decodeOneFrame(World world, final double[] codedSpectralEnvelope, final double[] frequencyAxis, int fftSize, final double[] melAxis, final FFT.FFTComplex weight, int maxDimension, int numberOfDimensions, final InverseComplexFFT inverseComplexFFT, double[] spectralEnvelope) {
        double[] melSpectrum = new double[maxDimension + 2];

        iDCTForCodec(world, codedSpectralEnvelope, maxDimension, weight, inverseComplexFFT, numberOfDimensions, melSpectrum, 1);
        melSpectrum[0] = melSpectrum[1];
        melSpectrum[maxDimension + 1] = melSpectrum[maxDimension];

        Utils.interp1(melAxis, melSpectrum, frequencyAxis, spectralEnvelope);

        for (int i = 0; i < fftSize / 2 + 1; ++i)
            spectralEnvelope[i] = Math.exp(spectralEnvelope[i] / maxDimension);
    }

    public static void getParametersForCoding(double floorFrequency, double ceilFrequency, int fs, int fftSize, double[] melAxis, double[] frequencyAxis, FFT.FFTComplex weight) {
        int maxDimension = fftSize / 2;
        double floorMel = frequencyToMel(floorFrequency);
        double ceilMel = frequencyToMel(ceilFrequency);

        for (int i = 0; i < maxDimension; ++i) {
            melAxis[i] = (ceilMel - floorMel) * i / maxDimension + floorMel;
            weight.get()[i][0] = 2.0 * Math.cos(i * Math.PI / fftSize) / Math.sqrt(fftSize);
            weight.get()[i][1] = 2.0 * Math.sin(i * Math.PI / fftSize) / Math.sqrt(fftSize);
        }
        weight.get()[0][0] /= Math.sqrt(2.0);

        for (int i = 0; i <= maxDimension; ++i)
            frequencyAxis[i] = frequencyToMel((double) (i) * fs / fftSize);
    }

    public static void getParametersForDecoding(double floorFrequency, double ceilFrequency, int fs, int fftSize, int numberOfDimensions, double[] melAxis, double[] frequencyAxis, FFT.FFTComplex weight) {
        int maxDimension = fftSize / 2;
        double floorMel = frequencyToMel(floorFrequency);
        double ceilMel = frequencyToMel(ceilFrequency);

        for (int i = 0; i < numberOfDimensions; ++i) {
            weight.get()[i][0] = Math.cos(i * Math.PI / fftSize) * Math.sqrt(fftSize);
            weight.get()[i][1] = Math.sin(i * Math.PI / fftSize) * Math.sqrt(fftSize);
        }
        weight.get()[0][0] /= Math.sqrt(2.0);
        for (int i = 0; i < maxDimension; ++i)
            melAxis[i + 1] =
                    melToFrequency((ceilMel - floorMel) * i / maxDimension + floorMel);
        melAxis[0] = 0;
        melAxis[maxDimension + 1] = fs / 2.0;

        for (int i = 0; i < fftSize / 2 + 1; ++i)
            frequencyAxis[i] = (double) (i) * fs / fftSize;
    }
}
