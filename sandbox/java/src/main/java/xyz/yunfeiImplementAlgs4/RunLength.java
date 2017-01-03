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
public class RunLength {   
    public static void compress(String filename, String output) {
        BinaryIn in = new BinaryIn(filename);
        BinaryOut out = new BinaryOut(output);
        int oneCount = 0;
        int zeroCount = 0;
        boolean firstZeroRunLength = false; //record whether first 8 bit for run length of 0 is written or not        
        while(!in.isEmpty()) {
            if (in.readBoolean()) {
                //make sure first 8 bit is for
                //run length of 0
                if (!firstZeroRunLength || zeroCount > 0) {
                    firstZeroRunLength = true;
                    out.write(zeroCount, 8);
                    zeroCount = 0;      
                }
                oneCount++;
                if(oneCount > 255) {
                    out.write(255, 8); //for 0
                    out.write(0, 8); //for 1
                    oneCount -= 255;                    
                }
            } else {
                if (oneCount > 0) {                    
                    out.write(oneCount, 8);
                    oneCount = 0;
                }
                zeroCount++;
                if (zeroCount > 255) {
                    out.write(255, 8); //for 1
                    out.write(0, 8);//for 0
                    zeroCount -= 255;
                }
            }
        }
        //make sure file ends with count of 1
        //such that we can read in pairs later
        if (oneCount > 0) {
            out.write(oneCount, 8);            
        } else {
            out.write(zeroCount, 8);
            out.write(0, 8);
        }
        out.close();
    }
    public static void expand(String filename, String output) {
        BinaryIn in = new BinaryIn(filename);
        BinaryOut out = new BinaryOut(output);
        int oneCount = 0;
        int zeroCount = 0;
        while(!in.isEmpty()) {
            //read 8 bits in pairs
            //first is for 0, the other for 1
            zeroCount = in.readByte();
            while(zeroCount-- > 0) {
                out.write(false);
            }
            oneCount = in.readByte();
            while(oneCount-- > 0) {
                out.write(true);
            }
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
