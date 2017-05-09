package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/15/17.
 */
/******************************************************************************
 *  Compilation:  javac KWIK.java
 *  Execution:    java KWIK file.txt
 *  Dependencies: StdIn.java StdOut.java In.java SuffixArray.java
 *  Data files:   http://algs4.cs.princeton.edu/63suffix/tale.txt
 *                http://algs4.cs.princeton.edu/63suffix/mobydick.txt
 *
 *  Keyword-in-context search.
 *
 *  %  java KWIK tale.txt 15
 *  majesty
 *   most gracious majesty king george th
 *  rnkeys and the majesty of the law fir
 *  on against the majesty of the people
 *  se them to his majestys chief secreta
 *  h lists of his majestys forces and of
 *
 *  the worst
 *  w the best and the worst are known to y
 *  f them give me the worst first there th
 *  for in case of the worst is a friend in
 *  e roomdoor and the worst is over then a
 *  pect mr darnay the worst its the wisest
 *  is his brother the worst of a bad race
 *  ss in them for the worst of health for
 *   you have seen the worst of her agitati
 *  cumwented into the worst of luck buuust
 *  n your brother the worst of the bad rac
 *   full share in the worst of the day pla
 *  mes to himself the worst of the strife
 *  f times it was the worst of times it wa
 *  ould hope that the worst was over well
 *  urage business the worst will be over i
 *  clesiastics of the worst world worldly
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.*;

/**
 *  The {@code KWIK} class provides a {@link SuffixArray} client for computing
 *  all occurrences of a keyword in a given string, with surrounding context.
 *  This is known as <em>keyword-in-context search</em>.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/63suffix">Section 6.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class KWIK {
  private KWIK() {}
  public static void main(String[] args) {
    In in = new In(args[0]);
    String q = args[1];
    String s = in.readAll().replaceAll("\\W+"," ");
    SuffixArrayInteface sa = null;
    sa = new SuffixArray(s);
    runBenchmark(sa, s, q);
    /*sa = new SuffixArrayX(s);
    runBenchmark(sa, s, q);*/
    /*sa = new Manber(s);
    runBenchmark(sa, s, q);*/
  }
  private static void runBenchmark(SuffixArrayInteface sa, String s, String q) {
    long totalTime = 0;
    for (int count = 0; count < 1; count++) {
      long startTime = System.nanoTime();
      int contextWidth = 15;
      int rank = sa.rank(q);
      if (rank == 0 || rank == s.length()) {
        System.out.println("cannot found");
      } else {
        System.out.println(s.substring(
                Math.max(0, sa.index(rank) - contextWidth),
                Math.min(s.length(), sa.index(rank) + q.length() + contextWidth)));
        for (int i = rank + 1; sa.lcp(i) >= q.length() && i < s.length(); i++) {
          System.out.println(s.substring(
                  Math.max(0, sa.index(i) - contextWidth),
                  Math.min(s.length(), sa.index(i) + q.length() + contextWidth)));
        }
      }
      long endTime = System.nanoTime();
      long duration = endTime - startTime;
      totalTime += duration;
      System.out.println("time elapsed: " + duration + " ms.");
    }
    System.out.println("total time: " + totalTime + " ms.");
  }
}
