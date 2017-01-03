package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;
/*
 * convert 000100111 to 3123, i.e. lengths of runs of bits
 *  in binary format (8 bits). Why not represent run length in
 *  more than 8 bits? because we don't expect such a long run of 
 *  1s or 0s (longer than 255).
 *  first 8 bit is run length for 0, then 1 and 0 alternates
 */
public class RunLength2 {   
    public static void compress(String filename, String output) {
        BinaryIn in = new BinaryIn(filename);
        BinaryOut out = new BinaryOut(output);
        int run = 0;
        boolean old = false;
        boolean current = false;
        while(!in.isEmpty()) {
            current = in.readBoolean();
            if (current != old) {
                out.write(run, 8);
                run = 1; //1 for the new bit
                old = !old;
            } else {
                if (run == 255) {
                    out.write(255, 8);
                    out.write(0, 8);
                    run = 0;
                }
                run++;
            }
        }
        out.write(run, 8);
        out.close();            
    }
    public static void expand(String filename, String output) {
        BinaryIn in = new BinaryIn(filename);
        BinaryOut out = new BinaryOut(output);
        int run = 0;
        boolean current = false;
        while(!in.isEmpty()) {
            //run = in.readByte(); //here we CAN'T use byte because it is signed!!! whereas run length is recorded by unsigned 8-bit int
            run = in.readInt(8);
            for (int i = 0; i < run; i++) {
                out.write(current);
            }
            current = !current;
        }
        out.close();
    }
    public static void main(String[] args) {
        //usage: 0 - filename outputname
        //usage: 0 + filename outputname
        if (args[0].equals("-")) {
            compress(args[1], args[2]);
        } else if (args[0].equals("+")) {
            expand(args[1], args[2]);
        } else {
            throw new IllegalArgumentException();
        }
    }

}
