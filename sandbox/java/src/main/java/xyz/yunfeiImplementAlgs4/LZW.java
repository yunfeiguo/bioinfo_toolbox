package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/21/17.
 */
/******************************************************************************
 *  Compilation:  javac LZW.java
 *  Execution:    java LZW - < input.txt   (compress)
 *  Execution:    java LZW + < input.txt   (expand)
 *  Dependencies: BinaryIn.java BinaryOut.java
 *  Data files:   http://algs4.cs.princeton.edu/55compression/abraLZW.txt
 *                http://algs4.cs.princeton.edu/55compression/ababLZW.txt
 *
 *  Compress or expand binary input from standard input using LZW.
 *
 *  WARNING: STARTING WITH ORACLE JAVA 6, UPDATE 7 the SUBSTRING
 *  METHOD TAKES TIME AND SPACE LINEAR IN THE SIZE OF THE EXTRACTED
 *  SUBSTRING (INSTEAD OF CONSTANT SPACE AND TIME AS IN EARLIER
 *  IMPLEMENTATIONS).
 *
 *  See <a href = "http://java-performance.info/changes-to-string-java-1-7-0_06/">this article</a>
 *  for more details.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;
import edu.princeton.cs.algs4.TST;

import java.util.HashMap;
import java.util.Map;

/**
 *  The {@code LZW} class provides static methods for compressing
 *  and expanding a binary input using LZW compression over the 8-bit extended
 *  ASCII alphabet with 12-bit codewords.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/55compress">Section 5.5</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class LZW {
  private static final int R = 4096; //# of elements in dictionary
  private static final int W = 12;  //width of compressed code
  private static TST<Integer> initializeTST(int n) {
    TST<Integer> st = new TST<>();
    for (int i = 0; i < n; i++) {
      st.put("" + (char) i, i);
    }
    return st;
  }
  public static void compress(BinaryIn in, BinaryOut out) {
    int count = 256;
    TST<Integer> st = initializeTST(count);
    String currentString = "";
    while(!in.isEmpty()) {
      currentString += in.readChar();
      if(st.contains(currentString)) {
        continue;
      } else {
        String prefix = st.longestPrefixOf(currentString);
        out.write(st.get(prefix), W);
        if(count < R - 1)
          st.put(currentString, count++);
        currentString = currentString.substring(prefix.length());
      }
    }
    out.write(st.get(currentString), W);
    out.write(R-1, W); //mark end
  }
  public static void expand(BinaryIn in, BinaryOut out) {
    int count = 256;
    String[] int2String = new String[R];
    for (int i = 0; i < count; i++) {
      int2String[i] = "" + (char) i;
    }
    String currentString = "";
    while(!in.isEmpty()) {
      int i = in.readInt(W);
      if (i == R - 1)
        break;
      if(int2String[i] == null) {
        //we see an unseen codeword, it must be the next
        //codeword that we'd like to add to the trie
        //therefore, the string must have same prefix as currentString
        //and the next string is just currentString+its 1st char
        int2String[count] = currentString + currentString.substring(0,1);
      }
      out.write(int2String[i]);
      if (count < R - 1 && currentString.length() > 0) {
        int2String[count++] = currentString + int2String[i].charAt(0);
      }
      currentString = int2String[i];
    }
  }
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
