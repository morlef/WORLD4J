package world.api;

public interface AudioInterface {
    void write(final double[] x, int fs, int nbit);

    int getAudioLength();

    double[] read();
}
