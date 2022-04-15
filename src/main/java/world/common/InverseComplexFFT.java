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
public class InverseComplexFFT {
    int fftSize;
    FFT.FFTComplex input;
    FFT.FFTComplex output;
    FFT.FFTPlan inverseFFT;
}
