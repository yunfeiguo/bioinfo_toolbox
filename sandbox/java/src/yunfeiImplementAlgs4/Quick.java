package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac Quick.java
 *  Execution:    java Quick < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/23quicksort/tiny.txt
 *                http://algs4.cs.princeton.edu/23quicksort/words3.txt
 *
 *  Sorts a sequence of strings from standard input using quicksort.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java Quick < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *       
 *  % java Quick < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *
 *
 *  Remark: For a type-safe version that uses static generics, see
 *
 *    http://algs4.cs.princeton.edu/23quicksort/QuickPedantic.java
 *
 ******************************************************************************/

import java.util.Random;

import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

/**
 *  The <tt>Quick</tt> class provides static methods for sorting an
 *  array and selecting the ith smallest element in an array using quicksort.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */

public class Quick {
    private Quick() {}
    public static void sort(Comparable[] a) {
        if (a == null) return;
        int n = a.length;        
        sort(a, 0, n - 1);
    }
    private static int partition(Comparable[] a, int lo, int hi) {    
        int pivot = (int) (Math.random() * (hi - lo + 1)) + lo;
        exch(a, pivot, lo);  //put pivot at the head of array
        int i = lo;
        int j = hi + 1;
        while (true) {
            //everything to the left of i (exclusive) should be no larger than pivot
            //everything to the right of j (exclusive) should be no smaller than pivot
            while (less(a[++i], a[lo])) //must break for a[lo] == a[++i], otherwise may run in quadratic time
                if (i == hi)
                    break;
            while (less(a[lo], a[--j])) ; //must break for a[lo] == a[--j], otherwise may run in quadratic time
            if (i >= j) break;
            else
                exch(a, i, j);
        }
       //why do we use j rather than i here?
        //because a[j]<=a[lo] is guaranteed, it is safe to switch a[j] and a[lo]
       exch(a, lo, j);   
        return j;
    }
    private static void sort(Comparable[] a, int lo, int hi) {
        if (lo >= hi) return;
        int pivot = partition(a, lo, hi);
        sort(a, lo, pivot - 1);
        sort(a, pivot + 1, hi);
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
         String[] a = {"O", "S"};
         //String[] a = {"A","O", "O","B"};
         //String[] a = {"S", "O", "R", "T", "E", "X", "A", "M", "P", "L", "E"};
         Quick.sort(a);
         show(a);
     }

}
