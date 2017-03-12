package yunfeiImplementAlgs4;

/******************************************************************************
 *  Compilation: javac MSD.java
 *  Execution:   java MSD < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/51radix/words3.txt
 *                http://algs4.cs.princeton.edu/51radix/shells.txt
 *
 *  Sort an array of strings or integers using MSD radix sort.
 *
 *  % java MSD < shells.txt
 *  are
 *  by
 *  sea
 *  seashells
 *  seashells
 *  sells
 *  sells
 *  she
 *  she
 *  shells
 *  shore
 *  surely
 *  the
 *  the
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code MSD} class provides static methods for sorting an
 *  array of extended ASCII strings or integers using MSD radix sort.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/51radix">Section 5.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class MSD {
  private static final int INTEGER_BITS = 32;
  private static final int BITS_PER_BYTE = 8;
  private static final int R = 256; //vocabulary size
  public static void sort(String[] a) {
    if (a == null) throw new IllegalArgumentException("input array null");
    sort(a, new String[a.length], 0, a.length - 1, 0);
  }

  /**
   * sort strings a[lo..hi] by i, i+1, ..., -th characters
   * @param a
   * @param lo
   * @param hi
   * @param i
   */
  private static void sort(String[] a, String[] aux, int lo, int hi, int i) {
    //to speed up for small arrays, we can use insertion sort
    if (lo > hi)
      return;
    int[] count = new int[R + 2];
    for (int k = lo; k <= hi; k++) {
      count[charAt(a[k], i) + 2]++;
    }
    for (int k = 1; k < count.length; k++) {
      count[k] += count[k - 1];
    }
    for (int k = lo; k <= hi; k++) {
      aux[lo + count[charAt(a[k], i) + 1]++] = a[k];
    }
    for (int k = lo; k <= hi; k++) {
      a[k] = aux[k];
    }
    /*
    sort strings with same prefixes by remaining characters
     */
    for (int k = 2; k < count.length; k++) {
      sort(a, aux, count[k - 1], count[k] - 1, i + 1);
    }
  }

  /**
   * @param s
   * @param i
   * @return ith char int if it exists, otherwise -1
   */
  private static int charAt(String s, int i) {
    if (i >= 0 && i < s.length())
      return s.charAt(i);
    return -1;
  }
  public static void sort(int[] a) {
    if (a == null) throw new IllegalArgumentException("input array null");
    int[] aux = new int[a.length];
    sort(a, aux, 0, a.length - 1, 0);
  }

  private static void sort(int[] a, int[] aux, int lo, int hi, int i) {
    if (lo > hi || i >= INTEGER_BITS/BITS_PER_BYTE)
      return;
    int[] count = new int[(int) Math.pow(2, BITS_PER_BYTE) + 1];
    for (int k = lo; k <= hi; k++) {
      count[getByte(a[k], i) + 1]++;
    }
    for (int k = 1; k < count.length; k++) {
      count[k] += count[k - 1];
    }
    for (int k = lo; k <= hi; k++) {
      aux[count[getByte(a[k], i)]++] = a[k];
    }
    for (int k = lo; k <= hi; k++) {
      a[k] = aux[k];
    }
    /*
    recursively sort each partition
     */
    sort(a, aux, 0, count[0] - 1, i + 1);
    for (int k = 1; k < count.length; k++) {
      sort(a, aux, count[k - 1], count[k] - 1, i+1);
    }
  }

  /**
   *
   * @param n
   * @param i 0~3
   * @return ith byte of 32-bit integer n
   */
  private static int getByte(int n, int i) {
    int mask = (1 << BITS_PER_BYTE) - 1;
    return (n >> (INTEGER_BITS - (i + 1)*BITS_PER_BYTE)) & mask;
  }
  /**
   * Reads in a sequence of extended ASCII strings from standard input;
   * MSD radix sorts them;
   * and prints them to standard output in ascending order.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    String[] a = in.readAllStrings();
    int n = a.length;
    sort(a);
    for (int i = 0; i < n; i++)
      StdOut.println(a[i]);
    int[] b = new int[]{8,1111,3,7,2222,6,5,0};
    sort(b);
    for (int i : b) {
      System.out.println(i);
    }
  }
}
