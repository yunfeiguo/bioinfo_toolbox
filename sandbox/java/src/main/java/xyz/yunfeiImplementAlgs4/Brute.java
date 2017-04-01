package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/1/17.
 */

import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac Brute.java
 *  Execution:    java Brute pattern text
 *  Dependencies: StdOut.java
 *
 *  Reads in two strings, the pattern and the input text, and
 *  searches for the pattern in the input text using brute force.
 *
 *  % java Brute abracadabra abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:               abracadabra
 *
 *  % java Brute rab abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:         rab
 *
 *  % java Brute rabrabracad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                        rabrabracad
 *
 *  % java Brute bcara abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                                   bcara
 *
 *  % java Brute abacad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern: abacad
 *
 ******************************************************************************/
public class Brute {
  /**
   * regular brute force
   * @param pattern
   * @param txt
   * @return
   */
  public static int search1(String pattern, String txt) {
    if (pattern.length() > txt.length()) throw new IllegalArgumentException("pattern is longer");
    //i is start of matching
    for (int i = 0; i < txt.length() - pattern.length() + 1; i++) {
      int j = 0;
      for (; j < pattern.length() && txt.charAt(i+j) == pattern.charAt(j); j++) { }
      if (j == pattern.length())
        return i;
    }
    return txt.length();
  }

  /**
   * back up version
   * @param pattern
   * @param txt
   * @return
   */
  public static int search2(String pattern, String txt) {
    if (pattern.length() > txt.length()) throw new IllegalArgumentException("pattern is longer");
    int j = 0;
    int i = 0;
    for (; i < txt.length() && j < pattern.length(); i++) {
      if (txt.charAt(i) == pattern.charAt(j)) {
        j++;
      } else {
        i -= j;
        j = 0;
      }
    }
    if (j == pattern.length())
      return i - j;
    return txt.length();
  }
  /****char array version*****/

  public static int search1(char[] pattern, char[] txt) {
    return search1(new String(pattern), new String(txt));
  }

  /**
   * back up version
   * @param pattern
   * @param txt
   * @return
   */
  public static int search2(char[] pattern, char[] txt) {
    return search2(new String(pattern), new String(txt));
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

    int offset1a = search1(pat, txt);
    int offset2a = search2(pat, txt);
    int offset1b = search1(pattern, text);
    int offset2b = search2(pattern, text);

    // print results
    StdOut.println("text:    " + txt);

    // from brute force search method 1a
    StdOut.print("pattern: ");
    for (int i = 0; i < offset1a; i++)
      StdOut.print(" ");
    StdOut.println(pat);

    // from brute force search method 2a
    StdOut.print("pattern: ");
    for (int i = 0; i < offset2a; i++)
      StdOut.print(" ");
    StdOut.println(pat);

    // from brute force search method 1b
    StdOut.print("pattern: ");
    for (int i = 0; i < offset1b; i++)
      StdOut.print(" ");
    StdOut.println(pat);

    // from brute force search method 2b
    StdOut.print("pattern: ");
    for (int i = 0; i < offset2b; i++)
      StdOut.print(" ");
    StdOut.println(pat);
  }
}
