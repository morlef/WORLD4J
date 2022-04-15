package world.common;

import lombok.AllArgsConstructor;
import lombok.Data;
import lombok.NoArgsConstructor;
import world.FFT;

@Data
@NoArgsConstructor
@AllArgsConstructor
public class InverseComplexFFT {
    int fftSize;
    FFT.FFTComplex input;
    FFT.FFTComplex output;
    FFT.FFTPlan inverseFFT;
}
