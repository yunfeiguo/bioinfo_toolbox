package yunfeiImplementAlgs4;

/******************************************************************************
 *  Compilation:  javac Quick3string.java
 *  Execution:    java Quick3string < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/51radix/words3.txt
 *                http://algs4.cs.princeton.edu/51radix/shells.txt
 *
 *  Reads string from standard input and 3-way string quicksort them.
 *
 *  % java Quick3string < shell.txt
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
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

/**
 *  The {@code Quick3string} class provides static methods for sorting an
 *  array of strings using 3-way radix quicksort.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/51radix">Section 5.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Quick3String {
  public static void sort(String[] a) {
    if (a == null) throw new NullPointerException("input string array null");
    StdRandom.shuffle(a);
    sort(a, 0, a.length - 1, 0);
  }

  /**
   * sort strings a[lo..hi] based on i, i+1, ..,-th strings
   * @param a
   * @param lo
   * @param hi
   * @param i
   */
  private static void sort(String[] a, int lo, int hi, int i) {
    if (lo >= hi)
      return;
    int pivot = lo;
    int small = lo + 1;
    int large = hi;
    int scan = lo + 1; //index of next unseen element
    /*
    rearrange the array to put smaller elements at the beginning
    equal elements in the middle, and larger elements at the end
    a[small..(scan - 1)]: equal to pivot
    a[lo..(small - 1)]: smaller than pivot
    a[(large + 1)...hi]: larger than pivot
     */
    while (scan <= large) {
      if (charAt(a[scan], i) < charAt(a[pivot], i)) {
        swap(a, small, scan);
        small++;
        scan++;
      } else if (charAt(a[scan], i) > charAt(a[pivot], i)) {
        swap(a, scan, large);
        large--;
      } else {
        scan++;
      }
    }
    swap(a, pivot, small - 1);
    small--;
    //recursively sort each partition
    sort(a, lo, small - 1, i);
    if (a[small].length() > i)
      //we only need to compare the equal elements (at ith position)
      //when there are unexamined chars
      sort(a, small, large, i + 1);
    sort(a, large + 1, hi, i);
  }
  private static void swap(Object[] a, int i, int j) {
    Object tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }

  /**
   * return ith char (in integer form) if it exists
   * otherwise -1
   * @param s
   * @param i
   * @return
   */
  private static int charAt(String s, int i) {
    if (i >= 0 && i < s.length())
      return s.charAt(i);
    return - 1;
  }
  /**
   * Reads in a sequence of fixed-length strings from standard input;
   * 3-way radix quicksorts them;
   * and prints them to standard output in ascending order.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {

    // read in the strings from standard input
    In in = new In(args[0]);
    String[] a = in.readAllStrings();
    int n = a.length;

    // sort the strings
    sort(a);

    // print the results
    for (int i = 0; i < n; i++)
      StdOut.println(a[i]);
  }
}
