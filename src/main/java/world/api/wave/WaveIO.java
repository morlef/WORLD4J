package world.api.wave;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.UnsupportedAudioFileException;

import world.api.AudioInterface;

public record WaveIO(File file) implements AudioInterface {
    @Override
    public void write(double[] x, int fs, int nbit) {
        byte[] buf = new byte[x.length * fs];

        for (int i = 0; i < x.length; i++) {
            ByteBuffer.wrap(buf, i * fs, fs).putDouble(x[i]);
        }

        try {
            AudioInputStream stream = new AudioInputStream(new ByteArrayInputStream(buf), new AudioFormat(AudioFormat.Encoding.PCM_FLOAT, fs, nbit, 1, 4, fs, false), (long) (fs * nbit));
            AudioSystem.write(stream, AudioFileFormat.Type.WAVE, file);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public int getAudioLength() {
        try {
            AudioInputStream stream = AudioSystem.getAudioInputStream(file);
            return (int) (stream.getFrameLength() * stream.getFormat().getFrameSize());
        } catch (UnsupportedAudioFileException | IOException e) {
            e.printStackTrace();
        }
        return 0;
    }

    @Override
    public double[] read() {
        try {
            AudioInputStream stream = AudioSystem.getAudioInputStream(file);
            byte[] array = new byte[(int) stream.getFrameLength() * stream.getFormat().getFrameSize()];
            stream.read(array);
            double[] buf = new double[(int) (stream.getFrameLength())];

            for (int i = 0; i < buf.length; i++) {
                buf[i] = ByteBuffer
                        .wrap(array, i * stream.getFormat().getFrameSize(), stream.getFormat().getFrameSize())
                        .getDouble();
            }
            return buf;
        } catch (UnsupportedAudioFileException | IOException e) {
            e.printStackTrace();
        }
        return new double[0];
    }

}
