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
public class ForwardRealFFT {
    int fftSize;
    double[] waveform;
    FFT.FFTComplex spectrum;
    FFT.FFTPlan forwardFFT;

    public ForwardRealFFT(int fftSize) {
        this.setFftSize(fftSize);
        this.setWaveform(new double[fftSize]);
        this.setSpectrum(new FFT.FFTComplex(fftSize));
        this.setForwardFFT(FFT.r2c(fftSize, this.getWaveform(), this.getSpectrum(), FFT.ESTIMATE));
    }
}
