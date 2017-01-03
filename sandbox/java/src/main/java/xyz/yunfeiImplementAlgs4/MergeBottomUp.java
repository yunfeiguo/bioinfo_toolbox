package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac MergeBU.java
 *  Execution:    java MergeBU < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/22mergesort/tiny.txt
 *                http://algs4.cs.princeton.edu/22mergesort/words3.txt
 *   
 *  Sorts a sequence of strings from standard input using
 *  bottom-up mergesort.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java MergeBU < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *    
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *  
 *  % java MergeBU < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>MergeBU</tt> class provides static methods for sorting an
 *  array using bottom-up mergesort.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class MergeBottomUp {
    
    //this class should not be instantiated
    private MergeBottomUp() {}
    public static void sort(Comparable[] a) {
        if (a == null || a.length == 1) return;
        //no recursion required
        int n = a.length;
        Comparable[] aux = new Comparable[n];
        for (int size = 1; size < n; size += size) {
            for (int i = 0; i < n - size; i += 2 * size) {
                //i < n - size: make sure at least one element in right subarray
                //i += 2*size: every time, we merge two subarrays of size size
                //Math.min(i + 2*size - 1, n - 1), make sure we do not go out of bound
                merge(a, aux, i, i + size - 1, Math.min(i + 2*size - 1, n - 1));
            }
        }
    }
    private static void merge(Comparable[] a, Comparable[] aux, int lo, int mid, int hi) {
        //merge is able to handle the case where two subarrays have very different lengths
        for (int i = lo; i <= hi; i++)
            aux[i] = a[i];
        int i = lo;
        int j = mid + 1;
        for (int k = lo; k <= hi; k++) {
            if (i > mid)                    a[k] = aux[j++];
            else if (j > hi)                a[k] = aux[i++];
            else if (less(aux[j], aux[i]))  a[k] = aux[j++];
            else                            a[k] = aux[i++]; //for stability
        }
    }
    /***********************************************************************
     *  Helper sorting functions.
     ***************************************************************************/
     
     // is v < w ?
     private static boolean less(Comparable v, Comparable w) {
         return v.compareTo(w) < 0;
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
         String[] a = {"S", "O", "R", "T", "E", "X", "A", "M", "P", "L", "E"};
         MergeBottomUp.sort(a);
         show(a);
     }

}
