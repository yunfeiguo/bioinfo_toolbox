package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac SuffixArray.java
 *  Execution:    java SuffixArray < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/63suffix/abra.txt
 *
 *  A data type that computes the suffix array of a string.
 *
 *   % java SuffixArray < abra.txt
 *    i ind lcp rnk  select
 *   ---------------------------
 *    0  11   -   0  "!"
 *    1  10   0   1  "A!"
 *    2   7   1   2  "ABRA!"
 *    3   0   4   3  "ABRACADABRA!"
 *    4   3   1   4  "ACADABRA!"
 *    5   5   1   5  "ADABRA!"
 *    6   8   0   6  "BRA!"
 *    7   1   3   7  "BRACADABRA!"
 *    8   4   0   8  "CADABRA!"
 *    9   6   0   9  "DABRA!"
 *   10   9   0  10  "RA!"
 *   11   2   2  11  "RACADABRA!"
 *
 *  See SuffixArrayX.java for an optimized version that uses 3-way
 *  radix quicksort and does not use the nested class Suffix.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

import java.util.Arrays;

/**
 *  The {@code SuffixArray} class represents a suffix array of a string of
 *  length <em>n</em>.
 *  It supports the <em>selecting</em> the <em>i</em>th smallest suffix,
 *  getting the <em>index</em> of the <em>i</em>th smallest suffix,
 *  computing the length of the <em>longest common prefix</em> between the
 *  <em>i</em>th smallest suffix and the <em>i</em>-1st smallest suffix,
 *  and determining the <em>rank</em> of a query string (which is the number
 *  of suffixes strictly less than the query string).
 *  <p>
 *  This implementation uses a nested class {@code Suffix} to represent
 *   a suffix of a string (using constant time and space) and
 *  {@code Arrays.sort()} to sort the array of suffixes.
 *  The <em>index</em> and <em>length</em> operations takes constant time
 *  in the worst case. The <em>lcp</em> operation takes time proportional to the
 *  length of the longest common prefix.
 *  The <em>select</em> operation takes time proportional
 *  to the length of the suffix and should be used primarily for debugging.
 *  <p>
 *  For alternate implementations of the same API, see
 *  {@link SuffixArrayX}, which is faster in practice (uses 3-way radix quicksort)
 *  and uses less memory (does not create {@code Suffix} objects)
 *  and <a href = "http://algs4.cs.princeton.edu/63suffix/SuffixArrayJava6.java.html">SuffixArrayJava6.java</a>,
 *  which relies on the constant-time substring extraction method that existed
 *  in Java 6.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/63suffix">Section 6.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */
public class SuffixArray {
  private Suffix[] suffixes;
  /**
   * Initializes a suffix array for the given {@code text} string.
   *
   * @param text the input string
   */
  public SuffixArray(String text) {
    if (text == null) throw new NullPointerException("input null");
    int n = text.length();
    suffixes = new Suffix[n];

    for (int i = 0; i < n; i++) {
      suffixes[i] = new Suffix(text, i);
    }
    Arrays.sort(suffixes);
  }

  private class Suffix implements Comparable<Suffix> {
    public String s;
    public int index;

    public Suffix(String s, int i) {
      if (i < 0 || i >= s.length()) throw new IndexOutOfBoundsException("");
      this.s = s;
      this.index = i;
    }

    public int length() {
      return s.length() - index;
    }

    public char charAt(int i) {
      return this.s.charAt(i + index);
    }

    @Override
    public int compareTo(Suffix that) {
      if (this == that)
        return 0;
      int n = Math.min(this.length(), that.length());
      for (int i = 0; i < n; i++) {
        int cmp = this.charAt(i) - that.charAt(i);
        if (cmp < 0) {
          return -1;
        } else if (cmp > 0) {
          return 1;
        }
      }
      return this.length() - that.length();
    }
    @Override
    public String toString() {
      return s.substring(index);
    }
  }

  /**
   * Returns the length of the input string.
   *
   * @return the length of the input string
   */
  public int length() {
    return suffixes.length;
  }

  /**
   * Returns the index into the original string of the <em>i</em>th smallest suffix.
   * That is, {@code text.substring(sa.index(i))} is the <em>i</em>th smallest suffix.
   *
   * @param i an integer between 0 and <em>n</em>-1
   * @return the index into the original string of the <em>i</em>th smallest suffix
   * @throws java.lang.IndexOutOfBoundsException unless {@code 0 <= i < n}
   */
  public int index(int i) {
    if (i < 0 || i >= length()) throw new IndexOutOfBoundsException("");
    return suffixes[i].index;
  }

  /**
   * Returns the length of the longest common prefix of the <em>i</em>th
   * smallest suffix and the <em>i</em>-1st smallest suffix.
   *
   * @param i an integer between 1 and <em>n</em>-1
   * @return the length of the longest common prefix of the <em>i</em>th
   * smallest suffix and the <em>i</em>-1st smallest suffix.
   * @throws java.lang.IndexOutOfBoundsException unless {@code 1 <= i < n}
   */
  public int lcp(int i) {
    if (i < 1 || i >= length()) throw new IndexOutOfBoundsException("no suffix before me");

    return lcp(suffixes[i - 1], suffixes[i]);
  }

  /**
   * compare 2 strings, and return longest common prefix length
   * @param a
   * @param b
   * @return
   */
  private int lcp(Suffix a, Suffix b) {
    int n = Math.min(a.length(), b.length());
    int i = 0;
    for (i = 0; i < n; i++) {
      if (a.charAt(i) != b.charAt(i))
        break;
    }
    return i;
  }

  /**
   * Returns the <em>i</em>th smallest suffix as a string.
   *
   * @param i the index
   * @return the <em>i</em> smallest suffix as a string
   * @throws java.lang.IndexOutOfBoundsException unless {@code 0 <= i < n}
   */
  public String select(int i) {
    if (i < 0 || i >= length()) throw new IndexOutOfBoundsException("");
    return suffixes[i].toString();
  }

  /**
   * Returns the number of suffixes strictly less than the {@code query} string.
   * We note that {@code rank(select(i))} equals {@code i} for each {@code i}
   * between 0 and <em>n</em>-1.
   *
   * @param query the query string
   * @return the number of suffixes strictly less than {@code query}
   */
  public int rank(String query) {
    if (query == null) throw new NullPointerException();
    int lo = 0;
    int hi = length() - 1;
    while (lo <= hi) {
      int mid = lo + (hi - lo)/2;
      int cmp = compare(query, suffixes[mid]);
      if (cmp < 0) {
        hi = mid - 1;
      } else if (cmp > 0) {
        lo = mid + 1;
      } else {
        return mid;
      }
    }
    return lo;
  }

  private int compare(String q, Suffix s) {
    int n = Math.min(q.length(), s.length());
    for (int i = 0; i < n; i++) {
      int cmp = q.charAt(i) - s.charAt(i);
      if (cmp < 0) {
        return -1;
      } else if (cmp > 0) {
        return 1;
      }
    }
    return q.length() - s.length();
  }

  /**
   * Unit tests the {@code SuffixArray} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    String s = in.readAll().replaceAll("\\W+", " ").trim();
    SuffixArray suffix = new SuffixArray(s);

    // StdOut.println("rank(" + args[0] + ") = " + suffix.rank(args[0]));

    StdOut.println("  i ind lcp rnk select");
    StdOut.println("---------------------------");

    for (int i = 0; i < s.length(); i++) {
      int index = suffix.index(i);
      String ith = "\"" + s.substring(index, Math.min(index + 50, s.length())) + "\"";
      assert s.substring(index).equals(suffix.select(i));
      int rank = suffix.rank(s.substring(index));
      if (i == 0) {
        StdOut.printf("%3d %3d %3s %3d %s\n", i, index, "-", rank, ith);
      } else {
        int lcp = suffix.lcp(i);
        StdOut.printf("%3d %3d %3d %3d %s\n", i, index, lcp, rank, ith);
      }
    }
  }
}
