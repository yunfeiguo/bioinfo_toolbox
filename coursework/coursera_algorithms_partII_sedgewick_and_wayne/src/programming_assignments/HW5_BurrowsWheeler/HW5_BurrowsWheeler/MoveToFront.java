package HW5_BurrowsWheeler;
import edu.princeton.cs.algs4.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
public class MoveToFront {
    private static final int R = 256; //alphabet size
    // apply move-to-front encoding, reading from standard input and writing to standard output
    public static void encode(String inputfile, String outputfile) {
        Logger logger = Logger.getLogger("MoveToFrontEncode");       
        long time = System.nanoTime();
        BinaryIn in = new BinaryIn(inputfile);
        BinaryOut out = new BinaryOut(outputfile);
        //use linkedlist to realize move to front operation
        LinkedList<Character> alphabet = new LinkedList<Character>();
        for (int i = 0; i < R; i++) {
            alphabet.add((char) i);
        }
        int count = 0;
        while(!in.isEmpty()) {
            char c = in.readChar();
            int i = alphabet.indexOf(c);
            out.write((char) i); //8-bit is sufficient for extended ASCII chars            
            if (i != 0) {
                alphabet.remove(i);
                alphabet.addFirst(c);
            }
            count++;
        }
        out.close();
        time = System.nanoTime() - time;
        logger.log(Level.INFO, "processed " + count + " ASCII chars in " + time + " nano seconds");
    }
    // apply move-to-front decoding, reading from standard input and writing to standard output
    public static void decode(String inputfile, String outputfile) {
        Logger logger = Logger.getLogger("MoveToFrontDecode");
        BinaryIn in = new BinaryIn(inputfile);
        BinaryOut out = new BinaryOut(outputfile);
        long time = System.nanoTime();
        int count = 0;
        LinkedList<Character> alphabet = new LinkedList<Character>();
        for (int i = 0; i < R; i++) {
            alphabet.add((char) i);
        }
        while(!in.isEmpty()) {
            int code = in.readChar();
            char c = alphabet.get(code);
            out.write(c);
            if (code > 0) {
                alphabet.remove(code);
                alphabet.addFirst(c);
            }
            count++;
        }
        out.close();
        time = System.nanoTime() - time;
        logger.log(Level.INFO, "processed " + count + " ASCII chars in " + time + " nano seconds");
    }

    // if args[0] is '-', apply move-to-front encoding
    // if args[0] is '+', apply move-to-front decoding
    public static void main(String[] args) {
        if (args[0].equals("-")) {
            encode(args[1], args[2]);
        } else if (args[0].equals("+")) {
            decode(args[1], args[2]);
        } else {
            throw new IllegalArgumentException();
        }
    }
}