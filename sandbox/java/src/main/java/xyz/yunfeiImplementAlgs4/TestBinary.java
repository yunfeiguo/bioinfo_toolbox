package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;
import edu.princeton.cs.algs4.BinaryStdOut;

public class TestBinary {
    public static void main(String[] args) {
        BinaryOut out = new BinaryOut(args[0]);
        out.write(false);
        out.close();
        /*BinaryIn in = new BinaryIn(args[0]);
        while(!in.isEmpty()) {
            if(in.readBoolean()) {
                System.out.print(1);
            } else {
                System.out.print(0);
            }
        }*/       
    }
}
