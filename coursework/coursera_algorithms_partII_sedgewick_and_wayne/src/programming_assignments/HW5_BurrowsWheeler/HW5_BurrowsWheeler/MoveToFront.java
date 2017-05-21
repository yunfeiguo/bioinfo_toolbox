package programming_assignments.HW5_BurrowsWheeler.HW5_BurrowsWheeler;

import edu.princeton.cs.algs4.*;

/**
 * Created by guoy28 on 5/7/17.
 */
public class MoveToFront {
  /** apply move-to-front encoding,
  * reading from standard input and writing to standard output
   **/
  private static int R = 256;
  private static int W = 8; //width
  public static void encode(BinaryIn in, BinaryOut out) {
    int[] coding = new int[R];
    int[] decoding = new int[R];
    for (int i = 0; i < R; i++) {
      coding[i] = i;
      decoding[i] = i;
    }
    while(!in.isEmpty()) {
      int c = in.readInt(W);
      out.write(coding[c], W);
      /*
      worst case: O(R*n)
      typical case: O(R+n)
       */
      for (int i = coding[c] - 1; i >= 0; i--) {
        coding[decoding[i]] = i + 1;
        decoding[i+1] = decoding[i];
      }
      coding[c] = 0;
      decoding[0] = c;
    }
  }

  // apply move-to-front decoding, reading from standard input and writing to standard output
  public static void decode(BinaryIn in, BinaryOut out) {
    int[] coding = new int[R];
    int[] decoding = new int[R];
    for (int i = 0; i < R; i++) {
      coding[i] = i;
      decoding[i] = i;
    }
    while(!in.isEmpty()) {
      int c = in.readInt(W);
      int d = decoding[c];
      out.write(decoding[c], W);
      for (int i = c - 1; i >= 0; i--) {
        //coding[decoding[i]] = i + 1;
        decoding[i + 1] = decoding[i];
      }
      //coding[d] = 0;
      decoding[0] = d;
    }
  }

  // if args[0] is '-', apply move-to-front encoding
  // if args[0] is '+', apply move-to-front decoding
  public static void main(String[] args) {
    BinaryIn in = new BinaryIn(args[1]);
    BinaryOut out = new BinaryOut(args[2]);
    if (args[0].compareTo("-") == 0) {
      encode(in, out);
    } else if (args[0].compareTo("+") == 0) {
      decode(in, out);
    } else {
      throw new IllegalArgumentException("only - or + allowed");
    }
    out.close();
  }
}
