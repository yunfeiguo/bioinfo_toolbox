package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/3/17.
 */
/******************************************************************************
 *  Compilation:  javac BoyerMoore.java
 *  Execution:    java BoyerMoore pattern text
 *  Dependencies: StdOut.java
 *
 *  Reads in two strings, the pattern and the input text, and
 *  searches for the pattern in the input text using the
 *  bad-character rule part of the Boyer-Moore algorithm.
 *  (does not implement the strong good suffix rule)
 *
 *  % java BoyerMoore abracadabra abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:               abracadabra
 *
 *  % java BoyerMoore rab abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:         rab
 *
 *  % java BoyerMoore bcara abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                                   bcara
 *
 *  % java BoyerMoore rabrabracad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                        rabrabracad
 *
 *  % java BoyerMoore abacad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern: abacad
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

import java.util.Arrays;

/**
 *  The {@code BoyerMoore} class finds the first occurrence of a pattern string
 *  in a text string.
 *  <p>
 *  This implementation uses the Boyer-Moore algorithm (with the bad-character
 *  rule, but not the strong good suffix rule).
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/53substring">Section 5.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */
public class BoyerMoore {
  private char[] pattern;
  private int[] rightmost; //rightmost index of each char

  public BoyerMoore(String pattern) {
    this(pattern.toCharArray(), 256);
  }
  /**
   * Preprocesses the pattern string.
   *
   * @param pattern the pattern string
   * @param R the alphabet size
   */
  public BoyerMoore(char[] pattern, int R) {
    /*
    construct skip heuristics
     */
    this.pattern = new char[pattern.length];
    //char[] is mutable!!!
    System.arraycopy(pattern, 0, this.pattern, 0, pattern.length);
    this.rightmost = new int[R];
    Arrays.fill(rightmost, -1);
    for (int i = 0; i < pattern.length; i++) {
      rightmost[pattern[i]] = i;
    }
  }
  public int search(String txt) {
    return search(txt.toCharArray());
  }
  public int search(char[] txt) {
    if (txt.length < pattern.length) throw new IllegalArgumentException("query string too short");
    for (int i = pattern.length - 1; i < txt.length;) {
      int j = pattern.length - 1;
      //compare from the rightmost position of pattern
      for (; j >= 0 && pattern[j] == txt[i - pattern.length + 1 + j]; j--) {}
      if (j < 0)
        return i - pattern.length + 1;
      else {
        //how much we advance if there is a mismatch
        /*
        aabcba
         caa
         ||
         \/
        aabcba
           caa
         */
        int skip = Math.max(1, j - rightmost[txt[i - pattern.length + 1 + j]]);
        i += skip;
      }
    }
    return txt.length;
  }
  /**
   * Takes a pattern string and an input string as command-line arguments;
   * searches for the pattern string in the text string; and prints
   * the first occurrence of the pattern string in the text string.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    String pat = args[0];
    String txt = args[1];
    char[] pattern = pat.toCharArray();
    char[] text    = txt.toCharArray();

    BoyerMoore boyermoore1 = new BoyerMoore(pat);
    BoyerMoore boyermoore2 = new BoyerMoore(pattern, 256);
    int offset1 = boyermoore1.search(txt);
    int offset2 = boyermoore2.search(text);

    // print results
    StdOut.println("text:    " + txt);

    StdOut.print("pattern: ");
    for (int i = 0; i < offset1; i++)
      StdOut.print(" ");
    StdOut.println(pat);

    StdOut.print("pattern: ");
    for (int i = 0; i < offset2; i++)
      StdOut.print(" ");
    StdOut.println(pat);
  }
}
