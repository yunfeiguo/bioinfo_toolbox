package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/5/17.
 */
/******************************************************************************
 *  Compilation:  javac KMP.java
 *  Execution:    java KMP pattern text
 *  Dependencies: StdOut.java
 *
 *  Reads in two strings, the pattern and the input text, and
 *  searches for the pattern in the input text using the
 *  KMP algorithm.
 *
 *  % java KMP abracadabra abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:               abracadabra
 *
 *  % java KMP rab abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:         rab
 *
 *  % java KMP bcara abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                                   bcara
 *
 *  % java KMP rabrabracad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                        rabrabracad
 *
 *  % java KMP abacad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern: abacad
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

import java.util.Arrays;

/**
 *  The {@code KMP} class finds the first occurrence of a pattern string
 *  in a text string.
 *  <p>
 *  This implementation uses a version of the Knuth-Morris-Pratt substring search
 *  algorithm. The version takes time as space proportional to
 *  <em>N</em> + <em>M R</em> in the worst case, where <em>N</em> is the length
 *  of the text string, <em>M</em> is the length of the pattern, and <em>R</em>
 *  is the alphabet size.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/53substring">Section 5.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */
public class KMP {
  /*
  DFA implementation
   */
  private int[][] dfa;
  private char[] pattern;
  private int l;

  public KMP(String pattern) {
    this(pattern.toCharArray(), 256);
  }
  public KMP(char[] pattern, int R) {
    this.pattern = pattern;
    this.l = pattern.length;
    this.dfa = new int[R][l];

    //fill DFA matrix
    //use slowIndex to track recurring sequence inside pattern
    int slowIndex = 0;
    for (char c = 0; c < R; c++) {
        dfa[c][0] = 0;
    }
    dfa[pattern[0]][0] = 1;

    for (int i = 1; i < l; i++) {
      for (char c = 0; c < R; c++) {
        dfa[c][i] = dfa[c][slowIndex];
      }
      dfa[pattern[i]][i] = i + 1;
      if (pattern[i] == pattern[slowIndex])
        slowIndex++;
      else
        slowIndex = 0;
    }
  }
  public int search(String txt) {
    return search(txt.toCharArray());
  }
  public int search(char[] txt) {
    if (txt.length < l) throw new IllegalArgumentException("query string too short");
    int matchIndex = 0;
    for (int i = 0; i < txt.length; i++) {
      matchIndex = dfa[txt[i]][matchIndex];
      if (matchIndex == l)
        return i - l + 1;
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

    KMP kmp1 = new KMP(pat);
    int offset1 = kmp1.search(txt);

    KMP kmp2 = new KMP(pattern, 256);
    int offset2 = kmp2.search(text);

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
