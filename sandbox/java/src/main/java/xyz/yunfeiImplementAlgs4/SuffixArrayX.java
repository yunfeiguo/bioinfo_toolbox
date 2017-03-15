package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/14/17.
 */
/******************************************************************************
 *  Compilation:  javac SuffixArrayX.java
 *  Execution:    java SuffixArrayX < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/63suffix/abra.txt
 *
 *  A data type that computes the suffix array of a string using 3-way
 *  radix quicksort.
 *
 *  % java SuffixArrayX < abra.txt
 *    i ind lcp rnk  select
 *  ---------------------------
 *    0  11   -   0  !
 *    1  10   0   1  A!
 *    2   7   1   2  ABRA!
 *    3   0   4   3  ABRACADABRA!
 *    4   3   1   4  ACADABRA!
 *    5   5   1   5  ADABRA!
 *    6   8   0   6  BRA!
 *    7   1   3   7  BRACADABRA!
 *    8   4   0   8  CADABRA!
 *    9   6   0   9  DABRA!
 *   10   9   0  10  RA!
 *   11   2   2  11  RACADABRA!
 *
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code SuffixArrayX} class represents a suffix array of a string of
 *  length <em>n</em>.
 *  It supports the <em>selecting</em> the <em>i</em>th smallest suffix,
 *  getting the <em>index</em> of the <em>i</em>th smallest suffix,
 *  computing the length of the <em>longest common prefix</em> between the
 *  <em>i</em>th smallest suffix and the <em>i</em>-1st smallest suffix,
 *  and determining the <em>rank</em> of a query string (which is the number
 *  of suffixes strictly less than the query string).
 *  <p>
 *  This implementation uses 3-way radix quicksort to sort the array of suffixes.
 *  For a simpler (but less efficient) implementations of the same API, see
 *  {@link SuffixArray}.
 *  The <em>index</em> and <em>length</em> operations takes constant time
 *  in the worst case. The <em>lcp</em> operation takes time proportional to the
 *  length of the longest common prefix.
 *  The <em>select</em> operation takes time proportional
 *  to the length of the suffix and should be used primarily for debugging.
 *  <p>
 *  This implementation uses '\0' as a sentinel and assumes that the charater
 *  '\0' does not appear in the text.
 *  <p>
 *  In practice, this algorithm runs very fast. However, in the worst-case
 *  it can be very poor (e.g., a string consisting of N copies of the same
 *  character. We do not shuffle the array of suffixes before sorting because
 *  shuffling is relatively expensive and a pathologial input for which
 *  the suffixes start out in a bad order (e.g., sorted) is likely to be
 *  a bad input for this algorithm with or without the shuffle.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/63suffix">Section 6.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */
public class SuffixArrayX {
  private static final int CUTOFF = 5;
  private int[] suffixIndices;
  private String s;
  /**
   * initialize {@code SuffixArrayX} object
   * @param s
   * @throws NullPointerException if input is null
   */
  public SuffixArrayX(String s) {
    if (s == null) throw new NullPointerException();
    this.s = s;
    this.suffixIndices = new int[s.length()];
    for (int i = 0; i < suffixIndices.length; i++) {
      suffixIndices[i] = i;
    }
    //no shuffle needed
    sort(s, suffixIndices, 0, suffixIndices.length - 1, 0);
  }

  /**
   * sort suffixes of s (indexed by indices[lo..hi]) based on
   * k,k+1,...,all the way to the last character
   * @param s
   * @param indices
   * @param lo
   * @param hi
   * @param k
   */
  private void sort(String s, int[] indices, int lo, int hi, int k) {
    if (lo >= hi)
      return;
    if (hi - lo + 1< CUTOFF) {
      insertion(s, indices, lo, hi, k);
      return;
    }
    int pivot = lo;
    int pivotChar = getChar(s, indices[pivot] + k);
    int small = lo + 1;
    int large = hi;
    int scan = small;

    while (scan <= large) {
      int currentChar = getChar(s, indices[scan] + k);
      if (currentChar > pivotChar) {
        swap(indices, scan, large--);
      } else if (currentChar < pivotChar) {
        swap(indices, scan++, small++);
      } else {
        scan++;
      }
    }
    swap(indices, --small, pivot);
    sort(s, indices, lo, small - 1, k);
    if (pivotChar >= 0)
      sort(s, indices, small, large, k + 1);
    sort(s, indices, large + 1, hi, k);
  }

  /**
   * use insertion sort to sort suffixes from
   * lo to hi based on k, k+1, k+2, ... characters
   */
  private void insertion(String s, int[] indices, int lo, int hi, int k) {
    for (int i = lo; i <= hi; i++) {
      for (int j = i + 1; j <= hi; j++) {
        int compareIndex = k;
        while(getChar(s, indices[i] + compareIndex) == getChar(s, indices[j] + compareIndex) &&
                indices[i] + compareIndex < s.length() &&
                indices[j] + compareIndex < s.length())
          compareIndex++;
        if (getChar(s, indices[i] + compareIndex) > getChar(s, indices[j] + compareIndex))
          swap(indices, i, j);
      }
    }
  }
  private void swap(int[] a, int i, int j) {
    int tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }

  /**
   *
   * @param s
   * @param i
   * @return ith char if exists, otherwise -1
   */
  private int getChar(String s, int i) {
    if (i >= 0 && i < s.length())
      return s.charAt(i);
    return -1;
  }

  /**
   *
   * @param i
   * @return original index of ith suffix
   */
  public int index(int i) {
    if (i < 0 || i >= suffixIndices.length) throw new IndexOutOfBoundsException();
    return suffixIndices[i];
  }

  /**
   *
   * @param q
   * @return number of suffixes exactly smaller than q
   */
  public int rank(String q) {
    int lo = 0;
    int hi = suffixIndices.length - 1;
    while (lo <= hi) {
      int mid = lo + (hi - lo)/2;
      int cmp = compare(s, mid, q);
      if (cmp < 0) {
        lo = mid + 1;
      } else if (cmp > 0) {
        hi = mid - 1;
      } else {
        return mid;
      }
    }
    return lo;
  }

  /**
   * compare s.substring(mid) with q
   * @param s
   * @param i
   * @param q
   */
  private int compare(String s, int i, String q) {
    int l = Math.min(s.length() - index(i), q.length());
    for (int j = 0; j < l; j++) {
      if (s.charAt(j + index(i)) > q.charAt(j)) {
        return 1;
      } else if (s.charAt(j + index(i)) < q.charAt(j)) {
        return -1;
      }
    }
    return s.length() - index(i) - q.length();
  }

  /**
   *
   * @param i
   * @return ith smallest suffix
   */
  public String select(int i) {
    if (i < 0 || i >= suffixIndices.length) throw new IndexOutOfBoundsException();
    return s.substring(suffixIndices[i]);
  }
  /**
   *
   * @param i
   * @return length of longest common prefix between i and i-1 suffixes
   */
  public int lcp(int i) {
    if (i < 1 || i >= suffixIndices.length) throw new IndexOutOfBoundsException();
    int l = Math.min(suffixIndices.length - index(i - 1), suffixIndices.length - index(i));
    int j = 0;
    for (; j < l; j++) {
      if (getChar(s, index(i - 1) + j) != getChar(s, index(i) + j)) {
        break;
      }
    }
    return j;
  }
  /**
   * Unit tests the {@code SuffixArrayx} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    String s = in.readAll().replaceAll("\\W+", " ").trim();
    SuffixArrayX suffix1 = new SuffixArrayX(s);
    SuffixArray suffix2 = new SuffixArray(s);
    boolean check = true;
    for (int i = 0; check && i < s.length(); i++) {
      if (suffix1.index(i) != suffix2.index(i)) {
        StdOut.println("suffix1(" + i + ") = " + suffix1.index(i));
        StdOut.println("suffix2(" + i + ") = " + suffix2.index(i));
        String ith = "\"" + s.substring(suffix1.index(i), Math.min(suffix1.index(i) + 50, s.length())) + "\"";
        String jth = "\"" + s.substring(suffix2.index(i), Math.min(suffix2.index(i) + 50, s.length())) + "\"";
        StdOut.println(ith);
        StdOut.println(jth);
        check = false;
      }
    }
    StdOut.println("  i ind lcp rnk  select");
    StdOut.println("---------------------------");

    for (int i = 0; i < s.length(); i++) {
      int index = suffix1.index(i);
      String ith = "\"" + s.substring(index, Math.min(index + 50, s.length())) + "\"";
      int rank = suffix1.rank(s.substring(index));
      assert s.substring(index).equals(suffix1.select(i));
      if (i == 0) {
        StdOut.printf("%3d %3d %3s %3d  %s\n", i, index, "-", rank, ith);
      }
      else {
        // int lcp  = suffix.lcp(suffix2.index(i), suffix2.index(i-1));
        int lcp  = suffix1.lcp(i);
        StdOut.printf("%3d %3d %3d %3d  %s\n", i, index, lcp, rank, ith);
      }
    }
  }
}
