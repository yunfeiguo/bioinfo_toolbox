package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/19/17.
 */
/******************************************************************************
 *  Compilation:  javac RunLength.java
 *  Execution:    java RunLength - < input.txt   (compress)
 *  Execution:    java RunLength + < input.txt   (expand)
 *  Dependencies: BinaryIn.java BinaryOut.java
 *  Data files:   http://algs4.cs.princeton.edu/55compression/4runs.bin
 *                http://algs4.cs.princeton.edu/55compression/q32x48.bin
 *                http://algs4.cs.princeton.edu/55compression/q64x96.bin
 *
 *  Compress or expand binary input from standard input using
 *  run-length encoding.
 *
 *  % java BinaryDump 40 < 4runs.bin
 *  0000000000000001111111000000011111111111
 *  40 bits
 *
 *  This has runs of 15 0s, 7 1s, 7 0s, and 11 1s.
 *
 *  % java RunLength - < 4runs.bin | java HexDump
 *  0f 07 07 0b
 *  4 bytes
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;

import static yunfeiImplementAlgs4.Genome.compress;
import static yunfeiImplementAlgs4.Genome.expand;

/**
 *  The {@code RunLength} class provides static methods for compressing
 *  and expanding a binary input using run-length coding with 8-bit
 *  run lengths.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/55compress">Section 5.5</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class RunLength {
  public static final int width = 8;
  public static final int maxCount = 255;
  private RunLength() {}
  public static void compress(BinaryIn in, BinaryOut out) {
    int count = 0;
    boolean currentB = false;
    while(!in.isEmpty()) {
      boolean b = in.readBoolean();
      if (b != currentB || count == maxCount) {
        out.write(count, width);
        if (b == currentB) {
          out.write(0, width);
        }
        currentB = b;
        count = 0;
      }
      count++;
    }
    out.write(count, width);
  }
  public static void expand(BinaryIn in, BinaryOut out) {
    boolean currentB = false;
    while(!in.isEmpty()) {
      int count = in.readInt(width);
      for (int i = 0; i < count; i++) {
        out.write(currentB);
      }
      currentB = !currentB;
    }
  }

  /**
   * Sample client that calls {@code compress()} if the command-line
   * argument is "-" an {@code expand()} if it is "+".
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    BinaryIn in = new BinaryIn(args[1]);
    BinaryOut out = new BinaryOut(args[2]);
    if      (args[0].equals("-")) {
      compress(in, out);
    } else if (args[0].equals("+")) {
      expand(in, out);
    } else throw new IllegalArgumentException("Illegal command line argument");
    out.close();
  }
}
