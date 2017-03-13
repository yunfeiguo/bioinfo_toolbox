package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/13/17.
 */
/******************************************************************************
 *  Compilation:  javac LongestRepeatedSubstring.java
 *  Execution:    java LongestRepeatedSubstring < file.txt
 *  Dependencies: StdIn.java SuffixArray.java
 *  Data files:   http://algs4.cs.princeton.edu/63suffix/tale.txt
 *                http://algs4.cs.princeton.edu/63suffix/tinyTale.txt
 *                http://algs4.cs.princeton.edu/63suffix/mobydick.txt
 *
 *  Reads a text string from stdin, replaces all consecutive blocks of
 *  whitespace with a single space, and then computes the longest
 *  repeated substring in that text using a suffix array.
 *
 *  % java LongestRepeatedSubstring < tinyTale.txt
 *  'st of times it was the '
 *
 *  % java LongestRepeatedSubstring < mobydick.txt
 *  ',- Such a funny, sporty, gamy, jesty, joky, hoky-poky lad, is the Ocean, oh! Th'
 *
 *  % java LongestRepeatedSubstring
 *  aaaaaaaaa
 *  'aaaaaaaa'
 *
 *  % java LongestRepeatedSubstring
 *  abcdefg
 *  ''
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code LongestRepeatedSubstring} class provides a {@link SuffixArray}
 *  client for computing the longest repeated substring of a string that
 *  appears at least twice. The repeated substrings may overlap (but must
 *  be distinct).
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/63suffix">Section 6.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  <p>
 *  See also {@link LongestCommonSubstring}.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class LongestRepeatedSubstring {
  private LongestRepeatedSubstring() {}

  public static String lrs(String s) {
    SuffixArray sa = new SuffixArray(s);
    int longest = -1;
    int longestIndex = -1;
    for (int i = 1; i < s.length(); i++) {
      int currentLengthCommonPrefix = sa.lcp(i);
      if (currentLengthCommonPrefix > longest) {
        longest = currentLengthCommonPrefix;
        longestIndex = i;
      }
    }
    return s.substring(sa.index(longestIndex), sa.index(longestIndex) + longest);
  }
  /**
   * Unit tests the {@code lrs()} method.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    String text = in.readAll().replaceAll("\\s+", " ");
    StdOut.println("'" + lrs(text) + "'");
  }
}
