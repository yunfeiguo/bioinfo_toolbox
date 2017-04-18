package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;

import java.util.HashMap;
import java.util.Map;

/******************************************************************************
 *  Compilation:  javac Alphabet.java
 *  Execution:    java Alphabet
 *  Dependencies: StdOut.java
 *  
 *  A data type for alphabets, for use with string-processing code
 *  that must convert between an alphabet of size R and the integers
 *  0 through R-1.
 *
 *  Warning: supports only the basic multilingual plane (BMP), i.e,
 *           Unicode characters between U+0000 and U+FFFF.
 *
 ******************************************************************************/

public class Alphabet {
  private Character[] int2Char;
  private Map<Character, Integer> char2Int;

  private static final Alphabet BINARY = new Alphabet("01");
  private static final Alphabet OCTAL = new Alphabet("01234567");
  private static final Alphabet DECIMAL = new Alphabet("0123456789");
  private static final Alphabet HEXADECIMAL = new Alphabet("0123456789ABCDEF");
  private static final Alphabet DNA = new Alphabet("ATCG");
  public static final Alphabet BASE64 = new Alphabet("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/");

  public Alphabet(String s) {
    int2Char = new Character[s.length()];
    char2Int = new HashMap<>();
    int charCount = 0;
    for (int i = 0; i < s.length(); i++) {
      int2Char[i] = s.charAt(i);
      char2Int.put(s.charAt(i), i);
    }
  }

    public int[] toIndices(String s) {
      int[] indices = new int[s.length()];
      for (int i = 0; i < s.length(); i++) {
        indices[i] = char2Int.get(s.charAt(i));
      }
      return indices;
    }
    public String toChars(int[] code) {
      StringBuilder sb = new StringBuilder();
      for (int i = 0; i < code.length; i++) {
        sb.append(int2Char[code[i]]);
      }
      return sb.toString();
    }

  /**
   * Unit tests the {@code Alphabet} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    int[]  encoded1 = Alphabet.BASE64.toIndices("NowIsTheTimeForAllGoodMen");
    String decoded1 = Alphabet.BASE64.toChars(encoded1);
    StdOut.println(decoded1);

    int[]  encoded2 = Alphabet.DNA.toIndices("AACGAACGGTTTACCCCG");
    String decoded2 = Alphabet.DNA.toChars(encoded2);
    StdOut.println(decoded2);

    int[]  encoded3 = Alphabet.DECIMAL.toIndices("01234567890123456789");
    String decoded3 = Alphabet.DECIMAL.toChars(encoded3);
    StdOut.println(decoded3);
  }
}