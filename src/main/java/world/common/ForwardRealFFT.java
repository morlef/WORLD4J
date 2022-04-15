package world.common;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import world.FFT;

@Data
@NoArgsConstructor
@AllArgsConstructor
public class ForwardRealFFT {
    int fftSize;
    double[] waveform;
    FFT.FFTComplex spectrum;
    FFT.FFTPlan forwardFFT;
}
