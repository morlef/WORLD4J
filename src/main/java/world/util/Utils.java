package world.util;

import world.FFT;
import world.World;
import world.common.ForwardRealFFT;
import world.common.InverseRealFFT;
import world.common.MinimumPhaseAnalysis;

import java.util.UUID;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/common.cpp, src/matlabfunctions.cpp
 */
public class Utils {
    public static double getSafeAperiodicity(double x) {
        return Math.max(0.001, Math.min(0.999999999999, x));
    }

    public static void filterForDecimate(final double[] x, int r, double[] y) {
        double[] a = new double[3], b = new double[2];
        switch (r) {
            case 11 -> {
                a[0] = 2.450743295230728;
                a[1] = -2.06794904601978;
                a[2] = 0.59574774438332101;
                b[0] = 0.0026822508007163792;
                b[1] = 0.0080467524021491377;
            }
            case 12 -> {
                a[0] = 2.4981398605924205;
                a[1] = -2.1368928194784025;
                a[2] = 0.62187513816221485;
                b[0] = 0.0021097275904709001;
                b[1] = 0.0063291827714127002;
            }
            case 10 -> {
                a[0] = 2.3936475118069387;
                a[1] = -1.9873904075111861;
                a[2] = 0.5658879979027055;
                b[0] = 0.0034818622251927556;
                b[1] = 0.010445586675578267;
            }
            case 9 -> {
                a[0] = 2.3236003491759578;
                a[1] = -1.8921545617463598;
                a[2] = 0.53148928133729068;
                b[0] = 0.0046331164041389372;
                b[1] = 0.013899349212416812;
            }
            case 8 -> {
                a[0] = 2.2357462340187593;
                a[1] = -1.7780899984041358;
                a[2] = 0.49152555365968692;
                b[0] = 0.0063522763407111993;
                b[1] = 0.019056829022133598;
            }
            case 7 -> {
                a[0] = 2.1225239019534703;
                a[1] = -1.6395144861046302;
                a[2] = 0.44469707800587366;
                b[0] = 0.0090366882681608418;
                b[1] = 0.027110064804482525;
            }
            case 6 -> {
                a[0] = 1.9715352749512141;
                a[1] = -1.4686795689225347;
                a[2] = 0.3893908434965701;
                b[0] = 0.013469181309343825;
                b[1] = 0.040407543928031475;
            }
            case 5 -> {
                a[0] = 1.7610939654280557;
                a[1] = -1.2554914843859768;
                a[2] = 0.3237186507788215;
                b[0] = 0.021334858522387423;
                b[1] = 0.06400457556716227;
            }
            case 4 -> {
                a[0] = 1.4499664446880227;
                a[1] = -0.98943497080950582;
                a[2] = 0.24578252340690215;
                b[0] = 0.036710750339322612;
                b[1] = 0.11013225101796784;
            }
            case 3 -> {
                a[0] = 0.95039378983237421;
                a[1] = -0.67429146741526791;
                a[2] = 0.15412211621346475;
                b[0] = 0.071221945171178636;
                b[1] = 0.21366583551353591;
            }
            case 2 -> {
                a[0] = 0.041156734567757189;
                a[1] = -0.42599112459189636;
                a[2] = 0.041037215479961225;
                b[0] = 0.16797464681802227;
                b[1] = 0.50392394045406674;
            }
            default -> {
                a[0] = 0.0;
                a[1] = 0.0;
                a[2] = 0.0;
                b[0] = 0.0;
                b[1] = 0.0;
            }
        }

        double[] w = {0.0, 0.0, 0.0};
        double wt;
        for (int i = 0; i < x.length; ++i) {
            wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];
            y[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];
            w[2] = w[1];
            w[1] = w[0];
            w[0] = wt;
        }
    }

    public static void fftShift(final double[] x, double[] y) {
        for (int i = 0; i < x.length / 2; ++i) {
            y[i] = x[i + x.length / 2];
            y[i + x.length / 2] = x[i];
        }
    }

    public static void histc(final double[] x, final double[] edges, int[] index) {
        int count = 1;

        int i = 0;
        for (; i < edges.length; ++i) {
            index[i] = 1;
            if (edges[i] >= x[0]) break;
        }
        for (; i < edges.length; ++i) {
            if (edges[i] < x[count]) {
                index[i] = count;
            } else {
                index[i--] = count++;
            }
            if (count == x.length) break;
        }
        count--;
        for (i++; i < edges.length; ++i) index[i] = count;
    }

    public static void interp1(final double[] x, final double[] y, double[] xi, double[] yi) {
        double[] h = new double[x.length - 1];
        int[] k = new int[xi.length];

        for (int i = 0; i < x.length - 1; ++i) h[i] = x[i + 1] - x[i];
        for (int i = 0; i < xi.length; ++i) {
            k[i] = 0;
        }

        histc(x, xi, k);

        for (int i = 0; i < x.length; ++i) {
            double s = (xi[i] - x[k[i] - 1]) / h[k[i] - 1];
            yi[i] = y[k[i] - 1] + s * (y[k[i]] - y[k[i] - 1]);
        }
    }

    public static void decimate(final double[] x, int r, double[] y) {
        final int kNFact = 9;
        double[] tmp1 = new double[x.length + kNFact * 2];
        double[] tmp2 = new double[x.length + kNFact * 2];

        for (int i = 0; i < kNFact; ++i) tmp1[i] = 2 * x[0] - x[kNFact - i];
        System.arraycopy(x, 0, tmp1, kNFact, kNFact + x.length - 9);
        for (int i = kNFact + x.length; i < 2 * kNFact + x.length; ++i)
            tmp1[i] = 2 * x[x.length - 1] - x[x.length - 2 - (i - (kNFact + x.length))];

        filterForDecimate(tmp1, r, tmp2);
        for (int i = 0; i < 2 * kNFact + x.length; ++i)
            tmp1[i] = tmp2[2 * kNFact + x.length - i - 1];
        filterForDecimate(tmp1, r, tmp2);
        for (int i = 0; i < 2 * kNFact + x.length; ++i)
            tmp1[i] = tmp2[2 * kNFact + x.length - i - 1];

        int nout = (x.length - 1) / r + 1;
        int nbeg = r - r * nout + x.length;

        int count = 0;
        for (int i = nbeg; i < x.length + kNFact; i += r)
            y[count++] = tmp1[i + kNFact - 1];
    }

    public static void diff(final double[] x, double[] y) {
        for (int i = 0; i < x.length - 1; ++i) y[i] = x[i + 1] - x[i];
    }

    public static void interp1Q(double x, double shift, final double[] y, final double[] xi, double[] yi) {
        double[] xiFraction = new double[xi.length];
        double[] deltaY = new double[y.length];
        int[] xiBase = new int[xi.length];

        for (int i = 0; i < xi.length; ++i) {
            xiBase[i] = (int) ((xi[i] - x) / shift);
            xiFraction[i] = (xi[i] - x) / shift - xiBase[i];
        }
        diff(y, deltaY);
        deltaY[y.length - 1] = 0.0;

        for (int i = 0; i < xi.length; ++i)
            yi[i] = y[xiBase[i]] + deltaY[xiBase[i]] * xiFraction[i];
    }

    public static void fastFFTFilt(World world, final double[] x, final double[] h, int fftSize, final ForwardRealFFT forwardRealFFT, final InverseRealFFT inverseRealFFT, double[] y) {
        FFT.FFTComplex xSpectrum = new FFT.FFTComplex();
        xSpectrum.set(new Double[fftSize][]);

        fastFFTFiltA(world, x, fftSize, forwardRealFFT);
        for (int i = 0; i <= fftSize / 2; ++i) {
            xSpectrum.get()[i][0] = forwardRealFFT.getSpectrum().get()[i][0];
            xSpectrum.get()[i][1] = forwardRealFFT.getSpectrum().get()[i][1];
        }

        fastFFTFiltA(world, h, fftSize, forwardRealFFT);

        for (int i = 0; i <= fftSize / 2; ++i) {
            inverseRealFFT.getSpectrum().get()[i][0] =
                    xSpectrum.get()[i][0] * forwardRealFFT.getSpectrum().get()[i][0] -
                            xSpectrum.get()[i][1] * forwardRealFFT.getSpectrum().get()[i][1];
            inverseRealFFT.getSpectrum().get()[i][1] =
                    xSpectrum.get()[i][0] * forwardRealFFT.getSpectrum().get()[i][1] +
                            xSpectrum.get()[i][1] * forwardRealFFT.getSpectrum().get()[i][0];
        }
        world.getFft().execute(inverseRealFFT.getInverseFFT());

        System.arraycopy(inverseRealFFT.getWaveform(), 0, y, 0, fftSize);
    }

    private static void fastFFTFiltA(World world, double[] x, int fftSize, ForwardRealFFT forwardRealFFT) {
        for (int i = 0; i < x.length; ++i)
            forwardRealFFT.getWaveform()[i] = x[i] / fftSize;
        for (int i = x.length; i < fftSize; ++i)
            forwardRealFFT.getWaveform()[i] = 0.0;
        world.getFft().execute(forwardRealFFT.getForwardFFT());
    }

    public static double std(final double[] x) {
        double average = 0.0;
        for (double v : x) average += v;
        average /= x.length;

        double s = 0.0;
        for (double v : x) s += Math.pow(v - average, 2.0);
        s /= (x.length - 1);

        return Math.sqrt(s);
    }

    public static void setParametersForLinearSmoothing(int boundary, int fftSize, int fs, double width, final double[] powerSpectrum, double[] mirroringSpectrum, double[] mirroringSegment, double[] frequencyAxis) {
        for (int i = 0; i < boundary; ++i)
            mirroringSpectrum[i] = powerSpectrum[boundary - i];
        if (fftSize / 2 + boundary - boundary >= 0)
            System.arraycopy(powerSpectrum, 0, mirroringSpectrum, boundary, fftSize / 2 + boundary - boundary);
        for (int i = fftSize / 2 + boundary; i <= fftSize / 2 + boundary * 2; ++i)
            mirroringSpectrum[i] = powerSpectrum[fftSize / 2 - (i - (fftSize / 2 + boundary))];

        mirroringSegment[0] = mirroringSpectrum[0] * fs / fftSize;
        for (int i = 1; i < fftSize / 2 + boundary * 2 + 1; ++i)
            mirroringSegment[i] = mirroringSpectrum[i] * fs / fftSize + mirroringSegment[i - 1];

        for (int i = 0; i <= fftSize / 2; ++i)
            frequencyAxis[i] = (double) i / fftSize * fs - width / 2.0;
    }

    public static int getSuitableFFTSize(int sample) {
        return (int) (Math.pow(2.0, (int) Math.log(sample) / World.LOG_2) + 1.0);
    }

    public static void dcCorrection(final double[] input, double f0, int fs, int fftSize, double[] output) {
        int upperLimit = 2 + (int) (f0 * fftSize / fs);
        double[] lowFrequencyReplica = new double[upperLimit];
        double[] lowFrequencyAxis = new double[upperLimit];

        for (int i = 0; i < upperLimit; ++i)
            lowFrequencyAxis[i] = (double) (i) * fs / fftSize;

        int upper_limit_replica = upperLimit - 1;
        interp1Q(f0 - lowFrequencyAxis[0], (double) (fs) / fftSize, input, lowFrequencyAxis, lowFrequencyReplica);

        for (int i = 0; i < upper_limit_replica; ++i)
            output[i] = input[i] + lowFrequencyReplica[i];
    }

    public static void linearSmoothing(final double[] input, double width, int fs, int fftSize, double[] output) {
        int boundary = (int) (width * fftSize / fs) + 1;

        double[] mirroringSpectrum = new double[fftSize / 2 + boundary * 2 + 1];
        double[] mirroringSegment = new double[fftSize / 2 + boundary * 2 + 1];
        double[] frequencyAxis = new double[fftSize / 2 + 1];
        setParametersForLinearSmoothing(boundary, fftSize, fs, width, input, mirroringSpectrum, mirroringSegment, frequencyAxis);

        double[] lowLevels = new double[fftSize / 2 + 1];
        double[] highLevels = new double[fftSize / 2 + 1];
        double originOfMirroringAxis = -(boundary - 0.5) * fs / fftSize;
        double discreteFrequencyInterval = (double) (fs) / fftSize;

        interp1Q(originOfMirroringAxis, discreteFrequencyInterval, mirroringSegment, frequencyAxis, lowLevels);

        for (int i = 0; i <= fftSize / 2; ++i) frequencyAxis[i] += width;

        interp1Q(originOfMirroringAxis, discreteFrequencyInterval, mirroringSegment, frequencyAxis, highLevels);

        for (int i = 0; i <= fftSize / 2; ++i)
            output[i] = (highLevels[i] - lowLevels[i]) / width;
    }

    public static void nuttallWindow(double[] y) {
        double tmp;
        for (int i = 0; i < y.length; ++i) {
            tmp  = i / (y.length - 1.0);
            y[i] = 0.355768 - 0.487396 * Math.cos(2.0 * Math.PI * tmp) + 0.144232 * Math.cos(4.0 * Math.PI * tmp) - 0.012604 * Math.cos(6.0 * Math.PI * tmp);
        }
    }

    public static void getMinimumPhaseSpectrum(World world, final MinimumPhaseAnalysis minimumPhase) {
        for (int i = minimumPhase.getFftSize() / 2 + 1; i < minimumPhase.getFftSize(); ++i)
            minimumPhase.getLogSpectrum()[i] = minimumPhase.getLogSpectrum()[minimumPhase.getFftSize() - i];

        world.getFft().execute(minimumPhase.getInverseFFT());
        minimumPhase.getCepstrum().get()[0][1] *= -1.0;
        for (int i = 1; i < minimumPhase.getFftSize() / 2; ++i) {
            minimumPhase.getCepstrum().get()[i][0] *= 2.0;
            minimumPhase.getCepstrum().get()[i][1] *= -2.0;
        }
        minimumPhase.getCepstrum().get()[minimumPhase.getFftSize() / 2][1] *= -1.0;
        for (int i = minimumPhase.getFftSize() / 2 + 1; i < minimumPhase.getFftSize(); ++i) {
            minimumPhase.getCepstrum().get()[i][0] = 0.0;
            minimumPhase.getCepstrum().get()[i][1] = 0.0;
        }

        world.getFft().execute(minimumPhase.getForwardFFT());

        double tmp;
        for (int i = 0; i <= minimumPhase.getFftSize() / 2; ++i) {
            tmp = Math.exp(minimumPhase.getMinimumPhaseSpectrum().get()[i][0] / minimumPhase.getFftSize());
            minimumPhase.getMinimumPhaseSpectrum().get()[i][0] = tmp * Math.cos(minimumPhase.getMinimumPhaseSpectrum().get()[i][1] / minimumPhase.getFftSize());
            minimumPhase.getMinimumPhaseSpectrum().get()[i][1] = tmp * Math.sin(minimumPhase.getMinimumPhaseSpectrum().get()[i][1] / minimumPhase.getFftSize());
        }
    }

    public static String generateRandomName() {
        return UUID.randomUUID().toString().substring(0, 8);
    }
}
