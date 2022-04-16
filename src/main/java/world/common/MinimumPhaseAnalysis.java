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
public class MinimumPhaseAnalysis {
    int fftSize;
    double[] logSpectrum;
    FFT.FFTComplex minimumPhaseSpectrum;
    FFT.FFTComplex cepstrum;
    FFT.FFTPlan inverseFFT;
    FFT.FFTPlan forwardFFT;

    public MinimumPhaseAnalysis(int fftSize) {
        this.setFftSize(fftSize);
        this.setLogSpectrum(new double[fftSize]);
        this.setMinimumPhaseSpectrum(new FFT.FFTComplex(fftSize));
        this.setCepstrum(new FFT.FFTComplex(fftSize));
        this.setInverseFFT(FFT.r2c(fftSize, this.getLogSpectrum(), this.getCepstrum(), FFT.ESTIMATE));
        this.setForwardFFT(FFT.dft(fftSize, this.getCepstrum(), this.getMinimumPhaseSpectrum(), FFT.FORWARD, FFT.ESTIMATE));
    }
}
