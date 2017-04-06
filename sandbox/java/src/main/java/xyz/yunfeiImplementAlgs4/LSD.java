package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac LSD.java
 *  Execution:    java LSD < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/51radix/words3.txt
 *
 *  LSD radix sort
 *
 *    - Sort a String[] array of n extended ASCII strings (R = 256), each of length w.
 *
 *    - Sort an int[] array of n 32-bit integers, treating each integer as
 *      a sequence of w = 4 bytes (R = 256).
 *
 *  Uses extra space proportional to n + R.
 *
 *
 *  % java LSD < words3.txt
 *  all
 *  bad
 *  bed
 *  bug
 *  dad
 *  ...
 *  yes
 *  yet
 *  zoo
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

import java.util.Arrays;

/**
 *  The {@code LSD} class provides static methods for sorting an
 *  array of <em>w</em>-character strings or 32-bit integers using LSD radix sort.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/51radix">Section 5.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class LSD {
  private static final int BITS_PER_BYTE = 8;
  private static final int INTEGER_BITS = 32;
  private LSD() {}

  /**
   * sort 32-bit integers using LSD
   * @param a
   */
  public static void sort(int[] a) {
    /*
    here if we treat each integer as 32-bit long binary string
    there is little speed gain

    if we partition each 32-bit binary string into 4-byte chars
    then we can be much faster
     */
    if (a == null) throw new NullPointerException("array null");
    int w = INTEGER_BITS / BITS_PER_BYTE;
    int R = (int) Math.pow(2, BITS_PER_BYTE); //vocabulary size
    int[] count = new int[R + 1];
    int[] aux = new int[a.length];

    for (int i = w - 1; i >= 0; i--) {
      Arrays.fill(count, 0);
      for (int j : a) {
        //convert ith byte to number (rightmost is 0)
        count[getByte(j, i) + 1]++;
      }
      for (int j = 1; j < count.length; j++) {
        count[j] = count[j] + count[j - 1];
      }
      for (int j : a) {
        aux[count[getByte(j, i)]++] = j;
      }
      for (int j = 0; j < a.length; j++) {
        a[j] = aux[j];
      }
    }
    return;
  }

  /**
   * return ith byte of integer n
   * @param n
   * @param i 0~3
   * @return
   */
  private static int getByte(int n, int i) {
    int mask = (1 << BITS_PER_BYTE) - 1;
    return ((n >> (INTEGER_BITS - 1 - i) * BITS_PER_BYTE) & mask);
  }

  /**
   * sort array of strings of length w
   * @param a string array
   * @param w # of chars per strings
   */
  public static void sort(String[] a, int w) {
    if (a == null) throw new NullPointerException("array null");
    int[] count = new int[(int) Math.pow(2, BITS_PER_BYTE) + 1];
    String[] aux = new String[a.length];

    for (int i = w - 1; i >= 0; i--) {
      Arrays.fill(count, 0);
      //count each char
      for (String s : a) {
        count[(int) s.charAt(i) + 1]++;
      }
      //cumsum for counts to determine boundaries
      for (int j = 1; j< count.length; j++) {
        count[j] = count[j] + count[j - 1];
      }
      //place each char
      for (String s : a) {
        aux[count[(int) s.charAt(i)]++] = s;
      }
      for (int k = 0; k < a.length; k++) {
        a[k] = aux[k];
      }
    }
  }
  /**
   * Reads in a sequence of fixed-length strings from standard input;
   * LSD radix sorts them;
   * and prints them to standard output in ascending order.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    String[] a = in.readAllStrings();
    int n = a.length;

    // check that strings have fixed length
    int w = a[0].length();
    for (int i = 0; i < n; i++)
      assert a[i].length() == w : "Strings must have fixed length";

    // sort the strings
    sort(a, w);

    // print results
    for (int i = 0; i < n; i++)
      StdOut.println(a[i]);

    int bN = 100000;
    int[] b = new int[bN];
    for (int i = 0; i < bN; i++) {
      b[i] = i;
    }
    StdRandom.setSeed(11L);
    StdRandom.shuffle(b);

    System.out.println("LSD sort");
    long startTime = System.nanoTime();
    LSD.sort(Arrays.copyOf(b, b.length));
    long endTime = System.nanoTime();
    System.out.println("elements sorted: " + b.length);
    System.out.println("Time used: " + (endTime - startTime) + " ms.");

    System.out.println("Arrays.sort()");
    startTime = System.nanoTime();
    Arrays.sort(Arrays.copyOf(b, b.length));
    endTime = System.nanoTime();
    System.out.println("elements sorted: " + b.length);
    System.out.println("Time used: " + (endTime - startTime) + " ms.");
  }
}
