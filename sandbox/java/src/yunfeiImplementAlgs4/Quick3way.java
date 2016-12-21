package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac Quick3way.java
 *  Execution:    java Quick3way < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/23quicksort/tiny.txt
 *                http://algs4.cs.princeton.edu/23quicksort/words3.txt
 *   
 *  Sorts a sequence of strings from standard input using 3-way quicksort.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java Quick3way < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *    
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *  
 *  % java Quick3way < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>Quick3way</tt> class provides static methods for sorting an
 *  array using quicksort with 3-way partitioning.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Quick3way {
    private Quick3way() {}
    public static void sort(Comparable[] a) {
        if (a == null) return;
        sort(a, 0, a.length - 1);
    }    
    //this is better when number of duplicates
    //is large, when it's smaller, use Bentley-Mcllroy 3-way partitioning
    private static void sort(Comparable[] a, int lo, int hi) {
        if (lo >= hi) return;
        //it is better to do partition and sorting in one place
        //we need 3 pointers here for small values, equal values and large values
        //respectively
        int pivot = (int) (Math.random() * (hi - lo + 1)) + lo;
        Comparable target = a[pivot]; //put pivot at the beginning
        int i = lo; //element to be examined
        int j = lo; //right boundary (exclusive) for small vals
        int k = hi; //left boundary (exclusive) for large vals
        while (i <= k) {
            int comp = target.compareTo(a[i]);
            if (comp > 0) exch(a, i++, j++);
            else if (comp < 0) exch(a, i, k--);
            else //equal
                i++;
        }
        sort(a, lo, j -1);
        sort(a, k + 1, hi);
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
         return isSorted(a, 0, a.length - 1);
     }

     private static boolean isSorted(Comparable[] a, int lo, int hi) {
         for (int i = lo + 1; i <= hi; i++)
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
        Quick3way.sort(a);
        show(a);
    }
}
