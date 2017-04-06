package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/5/17.
 */

import edu.princeton.cs.algs4.StdOut;

import java.util.Arrays;

/******************************************************************************
 *  Compilation:  javac KMPplus.java
 *  Execution:    java KMPplus pattern text
 *  Dependencies: StdOut.java
 *
 *  Knuth-Morris-Pratt algorithm over UNICODE alphabet.
 *
 *  % java KMPplus ABABAC BCBAABACAABABACAA
 *  text:    BCBAABACAABABACAA
 *  pattern:          ABABAC
 *
 *  % java KMPplus aabaaaba ccaabaabaabaaabaab
 *  text:    ccaabaabaabaaabaab
 *  pattern:         aabaaaba
 *
 *  % java KMPplus aabaaabb ccaabaabaabaaabaab
 *  text:    ccaabaabaabaaabaab
 *  pattern:                   aabaaabb
 *
 ******************************************************************************/
public class KMPplus {
  private int l;
  private int[] recurrentTracker;
  private String pattern;
  public KMPplus(String pattern) {
    this.l = pattern.length();
    this.pattern = pattern;
    recurrentTracker = new int[l];

    int recurrentIndex = 0;
    recurrentTracker[0] = 0;
    for (int i = 1; i < l; i++) {
      recurrentTracker[i] = recurrentIndex;
      if (pattern.charAt(i) == pattern.charAt(recurrentIndex)) {
        recurrentIndex++;
      } else
        recurrentIndex = 0;
    }
  }
  public int search(String text) {
    if (text.length() < l) throw new IllegalArgumentException("too short");
    int matchIndex = 0;
    for (int i = 0; i < text.length(); i++) {
      if (text.charAt(i) == pattern.charAt(matchIndex)) {
        matchIndex++;
      } else if (pattern.charAt(recurrentTracker[matchIndex]) == text.charAt(i)) {
        matchIndex = recurrentTracker[matchIndex] + 1;
      } else {
        matchIndex = 0;
      }
      if (matchIndex == l)
        return i - l + 1;
    }
    return text.length();
  }
  // test client
  public static void main(String[] args) {
    String pattern = args[0];
    String text    = args[1];
    int m = pattern.length();
    int n = text.length();

    // substring search
    KMPplus kmp = new KMPplus(pattern);
    int offset = kmp.search(text);

    // print results
    StdOut.println("m = " + m + ", n = " + n);
    StdOut.println("text:    " + text);
    StdOut.print("pattern: ");
    for (int i = 0; i < offset; i++)
      StdOut.print(" ");
    StdOut.println(pattern);
    StdOut.println("offset: " + offset);
  }
}
