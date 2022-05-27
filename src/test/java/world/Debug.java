package world;

import java.io.File;

public class Debug {
    public static void main(String[] args) {
        new Debug(new File(args[0]));
    }

    public Debug(File file) {

    }

    public void DisplayInformation(int fs, int nbit, int x_length) {
        World.getLogger().info("File information\n");
        World.getLogger().info("Sampling : %d Hz %d Bit\n", fs, nbit);
        World.getLogger().info("Length %d [sample]\n", x_length);
        World.getLogger().info("Length %f [sec]\n", (double) (x_length) / fs);
    }
}
