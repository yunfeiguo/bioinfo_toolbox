package programming_assignments.HW5_BurrowsWheeler.HW5_BurrowsWheeler2;
import java.io.File;
import java.io.IOException;
class BWCompress {
    public static void main(String[] args) throws IOException {
        //we need two temp files
        File mtf = File.createTempFile("tmp", "mtf");
        File bwt = File.createTempFile("tmp", "bwt");
        bwt.deleteOnExit();
        mtf.deleteOnExit();
        if (args[0].equals("-")) {        
            BurrowsWheeler.encode(args[1], bwt.getCanonicalFile().toString());
            MoveToFront.encode(bwt.getCanonicalFile().toString(), mtf.getCanonicalFile().toString());
            Huffman.compress(mtf.getCanonicalPath().toString(), args[2]);
        } else if (args[0].equals("+")) {
            Huffman.expand(args[1], mtf.getCanonicalFile().toString());
            MoveToFront.decode(mtf.getCanonicalFile().toString(), bwt.getCanonicalFile().toString());
            BurrowsWheeler.decode(bwt.getCanonicalFile().toString(), args[2]);
        } else
            throw new UnsupportedOperationException();
    }
}
