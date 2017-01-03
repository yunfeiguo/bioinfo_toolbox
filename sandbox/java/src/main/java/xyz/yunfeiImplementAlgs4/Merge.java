/******************************************************************************
 *  Compilation:  javac Merge.java
 *  Execution:    java Merge < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/22mergesort/tiny.txt
 *                http://algs4.cs.princeton.edu/22mergesort/words3.txt
 *   
 *  Sorts a sequence of strings from standard input using mergesort.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java Merge < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *    
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *  
 *  % java Merge < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *  
 ******************************************************************************/
package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>Merge</tt> class provides static methods for sorting an
 *  array using mergesort.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/22mergesort">Section 2.2</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  For an optimized version, see {@link MergeX}.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Merge {
    private Merge() {} //this class should not be instantiated
    //stably merge a[lo..mid] with a[mid+1..hi] using aux[lo..hi]
    public static void merge(Comparable[] a, Comparable[] aux, int lo, int mid, int hi) {
        //precondition
        assert isSorted(a, lo, mid);
        assert isSorted(a, mid + 1, hi);
        if (lo >= hi) return;
        //copy to aux[]
        for (int i = lo; i <= hi; i++) {
            aux[i] = a[i];
        }
        //merge two halves
        int i = lo;
        int j = mid + 1;
        for (int k = lo; k <= hi; k++) {
            if      (j > hi)    a[k] = aux[i++];
            else if (i > mid)   a[k] = aux[j++];
            else if (less(aux[i], aux[j])) a[k] = aux[i++];
            else                       a[k] = aux[j++];
        }
        //postcondition
        assert isSorted(a, lo, hi);
    }
    public static void sort(Comparable[] a) {
        Comparable[] aux = new Comparable[a.length];
        sort(a, aux, 0, a.length - 1); //overload
    }
    public static void sort(Comparable[] a, Comparable[] aux, int lo, int hi) {
        //precondition
        assert a.length == aux.length;
        if (hi <= lo) return;
        int mid = lo + (hi - lo) / 2;
        sort(a, aux, lo, mid);
        sort(a, aux, mid + 1, hi);
        merge(a, aux, lo, mid, hi);
        //postcondition
        assert isSorted(a, lo, hi);
    }
    private static boolean less(Comparable x, Comparable y) {
        return x.compareTo(y) < 0;
    }
    private static boolean isSorted(Comparable[] a, int lo, int hi) {
        for (int i = lo + 1; i <= hi; i++)
            if (less(a[i], a[i - 1])) return false;
        return true;
    }
    private static void show(Comparable[] a) {
        for (int i = 0; i < a.length; i++) {
            StdOut.println(a[i]);
        }
    }
    public static void main(String[] args) {
        //String[] a = {};
        String[] a = {"A"};
        //String[] a = {"S", "O", "R", "T", "E", "X", "A", "M", "P", "L", "E"};
        Merge.sort(a);
        show(a);
    }

}
