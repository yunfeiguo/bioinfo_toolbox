/******************************************************************************
 *  Compilation:  javac QuickX.java
 *  Execution:    java QuickX N
 *  Dependencies: StdOut.java StdIn.java
 *  
 *  Uses the Bentley-McIlroy 3-way partitioning scheme,
 *  chooses the partitioning element using Tukey's ninther,
 *  and cuts off to insertion sort.
 *
 *  Reference: Engineering a Sort Function by Jon L. Bentley
 *  and M. Douglas McIlroy. Softwae-Practice and Experience,
 *  Vol. 23 (11), 1249-1265 (November 1993).
 *
 ******************************************************************************/

/**
 *  The <tt>QuickX</tt> class provides static methods for sorting an
 *  array using an optimized version of quicksort (using Bentley-McIlroy
 *  3-way partitioning, Tukey's ninther, and cutoff to insertion sort).
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;

public class QuickX {
    private static final int CUTOFF = 7; //cutoff for insertion sort
    private QuickX() {}
    public static void sort(Comparable[] a) {
        if (a == null) return;
        sort(a, 0, a.length - 1);
    }
    private static void sort(Comparable[] a, int lo, int hi) {
        if (lo >= hi) return;
        if (hi - lo <= CUTOFF) insertionSort(a, lo, hi);
        /* we need four pointers
         * |====|<<<<<<|???????????|>>>>>|=====|
         * lo   k      i->       <-j     l     hi
         * k: right boundary (exclusive) for equal values
         * i: element to be examined for left half
         * j: element to be examined for right half
         * l: left boundary (exclusive) for equal values
         * [k,i), values smaller than target
         * (j, l], values larger than target
         *
         **/
        Comparable target = a[(int) (Math.random()*(hi - lo + 1)) + lo];
        int k = lo;
        int i = lo - 1;
        int l = hi;
        int j = hi + 1;
        while (true) {
            while (less(a[++i], target))
                if (i == j)
                    break;
            while (less(target, a[--j]))
                if (i == j)
                    break;
            if (i == j && eq(a[i], target))
                exch(a, k++, i); //put last equal item to the beginning
            if (i >= j) break;
            
            exch(a, i, j); //doesn't matter if it's a[i, j] == target or a[i, j] < target
            if (eq(a[i], target))
                exch(a, i, k++);
            if (eq(a[j], target))
                exch(a, j, l--);
        }
        //then we move equal values to the center
        if (i == j) j++;
        while (k > lo)
            exch(a, --k, i--);        
        while (l < hi)
            exch(a, ++l, j++);        
        sort(a, lo, i);
        sort(a, j, hi);    
    }
    // sort from a[lo] to a[hi] using insertion sort
    private static void insertionSort(Comparable[] a, int lo, int hi) {
        for (int i = lo; i <= hi; i++)
            for (int j = i; j > lo && less(a[j], a[j-1]); j--)
                exch(a, j, j-1);
    }


    // return the index of the median element among a[i], a[j], and a[k]
    private static int median3(Comparable[] a, int i, int j, int k) {
        return (less(a[i], a[j]) ?
               (less(a[j], a[k]) ? j : less(a[i], a[k]) ? k : i) :
               (less(a[k], a[j]) ? j : less(a[k], a[i]) ? k : i));
    }

   /***************************************************************************
    *  Helper sorting functions.
    ***************************************************************************/
    
    // is v < w ?
    private static boolean less(Comparable v, Comparable w) {
        return v.compareTo(w) < 0;
    }

    // does v == w ?
    private static boolean eq(Comparable v, Comparable w) {
        return v.compareTo(w) == 0;
    }
        
    // exchange a[i] and a[j]
    private static void exch(Object[] a, int i, int j) {
        Object swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }


   /***************************************************************************
    *  Check if array is sorted - useful for debugging.
    ***************************************************************************/
    private static boolean isSorted(Comparable[] a) {
        for (int i = 1; i < a.length; i++)
            if (less(a[i], a[i-1])) return false;
        return true;
    }

    // print array to standard output
    private static void show(Comparable[] a) {
        for (int i = 0; i < a.length; i++) {
            StdOut.println(a[i]);
        }
    }
    public static void main(String[] args) {
        //String[] a = {"S", "O"};
        //String[] a = {"O", "S"};
        //String[] a = {"A","O", "O","B"};
        String[] a = {"O", "S", "O", "R", "T", "E", "E", "E", "X", "A", "M", "E", "P", "L", "E"};
        Quick.sort(a);
        show(a);
    }
}
