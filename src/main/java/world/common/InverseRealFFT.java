package world.common;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import world.FFT;

/**
 * Copyright 2012 Masanori Morise
 * @author mmorise [at] meiji.ac.jp (Masanori Morise)
 *
 * from src/world/common.h
 */
@Data
@NoArgsConstructor
@AllArgsConstructor
public class InverseRealFFT {
    int fftSize;
    double[] waveform;
    FFT.FFTComplex spectrum;
    FFT.FFTPlan inverseFFT;

    public InverseRealFFT(int fftSize) {
        this.setFftSize(fftSize);
        this.setWaveform(new double[fftSize]);
        this.setSpectrum(new FFT.FFTComplex(fftSize));
        this.setInverseFFT(FFT.c2r(fftSize, this.getSpectrum(), this.getWaveform(), FFT.ESTIMATE));
    }
}
