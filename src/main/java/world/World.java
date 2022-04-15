package world;

import lombok.Getter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import world.synthesis.RealtimeSynthesis;
import world.synthesis.Synthesis;
import world.synthesis.WorldSynthesizer;
import world.util.Utils;

public class World extends Thread {
    @Getter
    private static final Logger logger;

    public static final double CUT_OFF = 50.0;

    public static final double FLOOR_F0_STONEMASK = 40.0;

    public static final double MY_SAFE_GUARD_MINIMUM = 0.000000000001;
    public static final double EPS = 0.00000000000000022204460492503131;
    public static final double FLOOR_F0 = 71.0;
    public static final double CEIL_F0 = 800.0;
    public static final double DEFAULT_F0 = 500.0;
    public static final double LOG_2 = 0.69314718055994529;

    public static final double MAXIMUM_VALUE = 100000.0;

    public static final int HANNING = 1;
    public static final int BLACKMAN = 2;
    public static final double FREQUENCY_INTERVAL = 3000.0;
    public static final double UPPER_LIMIT = 15000.0;
    public static final double THRESHOLD = 0.85;
    public static final double FLOOR_F0D4C = 47.0;

    public static final double M0 = 1127.01048;
    public static final double F0 = 700.0;
    public static final double FLOOR_FREQUENCY = 40.0;
    public static final double CEIL_FREQUENCY = 20000.0;

    @Getter private final WorldSynthesizer synthesizer;

    @Getter private final Synthesis synthesis;
    @Getter private final RealtimeSynthesis realtimeSynthesis;

    @Getter private final CheapTrick cheapTrick;
    @Getter private final Codec codec;
    @Getter private final D4C d4c;
    @Getter private final Dio dio;
    @Getter private final FFT fft;
    @Getter private final Harvest harvest;
    @Getter private final StoneMask stoneMask;


    public World(WorldSynthesizer synthesizer) {
        super("WORLD/" + Utils.generateRandomName());

        this.synthesizer = synthesizer;

        synthesis = new Synthesis(this);
        realtimeSynthesis = new RealtimeSynthesis(this);

        cheapTrick = new CheapTrick(this);
        codec = new Codec(this);
        d4c = new D4C(this);
        dio = new Dio(this);
        fft = new FFT();
        harvest = new Harvest(this);
        stoneMask = new StoneMask(this);
    }

    @Override
    public void start() {
        super.start();
    }

    static {
        logger = LogManager.getLogger("WORLD");
    }
}
