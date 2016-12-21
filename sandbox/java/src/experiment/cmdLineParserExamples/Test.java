package experiment.cmdLineParserExamples;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
/**
 * Created by guoy28 on 9/22/16.
 */
public class Test {

    @Option(name = "-s" )
    public Stock s;

    /**
     * try run it with command line arguments "-s aapl" and "-s xx"
     * second should fail
     * @param args
     */
    public static void main(String[] args) {
        Test t = new Test();
        CmdLineParser parser = new CmdLineParser(t);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.out.println(e.getMessage());
        }
        System.out.println(t.s);
    }
}
