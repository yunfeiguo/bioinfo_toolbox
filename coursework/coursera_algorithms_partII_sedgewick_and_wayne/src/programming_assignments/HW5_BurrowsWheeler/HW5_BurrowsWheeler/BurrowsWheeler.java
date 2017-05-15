package programming_assignments.HW5_BurrowsWheeler.HW5_BurrowsWheeler;

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.Out;

/**
 * Created by guoy28 on 5/7/17.
 */
public class BurrowsWheeler {
  private static final int R = 256; //alphabet size
  // apply Burrows-Wheeler encoding, reading from standard input and writing to standard output
  public static void encode(BinaryIn in, BinaryOut out) {
    String s = in.readString();
    CircularSuffixArray sa = new CircularSuffixArray(s);
    for (int i = 0; i < s.length(); i++) {
      if (sa.index(i) == 0) {
        out.write(i);
      }
    }
    for (int i = 0; i < s.length(); i++) {
      out.write(s.charAt(sa.index(i)));
    }
  }

  // apply Burrows-Wheeler decoding, reading from standard input and writing to standard output
  public static void decode(BinaryIn in, BinaryOut out) {
    int first = in.readInt();
    //if we sort all circular suffix arrays
    //the encoded msg is the last column
    char[] lastColumn = in.readString().toCharArray();
    int n = lastColumn.length;
    int[] next = new int[n];
    //sort by key-indexed counting
    //goal is stable sorting
    int[] counts = new int[R + 1];
    //count occurrences
    for (char c : lastColumn) {
      counts[c + 1]++;
    }
    //calc boundary
    for (int i = 1; i < counts.length; i++) {
      counts[i] += counts[i - 1];
    }
    for (int i = 0; i < n; i++) {
      next[counts[lastColumn[i]]++] = i;
    }
    //restore original string
    for (int i = 0; i < n; i++) {
      out.write(lastColumn[next[first]]);
      first = next[first];
    }
  }

  // if args[0] is '-', apply Burrows-Wheeler encoding
  // if args[0] is '+', apply Burrows-Wheeler decoding
  public static void main(String[] args) {
    BinaryIn in = new BinaryIn(args[1]);
    BinaryOut out = new BinaryOut(args[2]);
    if (args[0].compareTo("-") == 0) {
      decode(in, out);
    } else if (args[0].compareTo("+") == 0) {
      encode(in, out);
    } else {
      throw new IllegalArgumentException("only - or + allowed");
    }
    out.close();
  }
}
