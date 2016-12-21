/******************************************************************************
 *  Compilation:  javac MergeX.java
 *  Execution:    java MergeX < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/22mergesort/tiny.txt
 *                http://algs4.cs.princeton.edu/22mergesort/words3.txt
 *   
 *  Sorts a sequence of strings from standard input using an
 *  optimized version of mergesort.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java MergeX < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *    
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *  
 *  % java MergeX < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Stopwatch;
/**
 *  The <tt>MergeX</tt> class provides static methods for sorting an
 *  array using an optimized version of mergesort.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/22mergesort">Section 2.2</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class MergeX {
    //speed up mergesort with insertion sort
    private static final int CUTOFF = 7; //cutoff for insertion sort
    private MergeX() {}
    public static void sort(Comparable[] a) {
        if (a == null) return;
        Comparable[] aux = new Comparable[a.length];
        sort(a, aux, 0, a.length - 1);
    }
    public static void sort(Comparable[] a, Comparable[] aux, int lo, int hi) {
        if (lo >= hi) return;
        if (hi - lo <= CUTOFF) {
            insertion(a, lo, hi);
            return;
        }
        int mid = lo + (hi - lo) / 2;
        sort(a, aux, lo, mid);
        sort(a, aux, mid + 1, hi);
        merge(a, aux, lo, mid, hi);
    }
    private static void merge(Comparable[] a, Comparable[] aux, int lo, int mid, int hi) {
        //copy elements to auxillary array
        for (int i = lo; i <= hi; i++)
            aux[i] = a[i];
        int i = lo; //pointer for first half
        int j = mid + 1; //pointer for second half
        for (int k = lo; k <= hi; k++) {
            if (i > mid)        a[k] = aux[j++];
            else if (j > hi)    a[k] = aux[i++];
            else if (less(aux[j], aux[i])) a[k] = aux[j++];
            //this is critical for stable sorting
            //if we use aux[j++], original order will be broken
            else                       a[k] = aux[i++];
        }        
    }
    private static void insertion(Comparable[] a, int lo, int hi) {
        if (lo >= hi) return;
        for (int i = lo; i <= hi; i++) {
            for (int j = i - 1; j >= lo && less(a[j + 1], a[j]); j--) {
                exch(a, j + 1, j);
            }
        }
    }
           
    // exchange a[i] and a[j]
    private static void exch(Object[] a, int i, int j) {
        Object swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }

    // is a[i] < a[j]?
    private static boolean less(Comparable a, Comparable b) {
        return a.compareTo(b) < 0;
    }
    private static void show(Object[] a) {
        for (int i = 0; i < a.length; i++) {
            StdOut.println(a[i]);
        }
    }   
    public static void main(String[] args) {
        String[] a = StdIn.readAllStrings();
        String[] b = new String[a.length];
        for (int i = 0; i < a.length; i++) {
            b[i] = a[i];
        }
        Stopwatch timer = new Stopwatch();
        //edu.princeton.cs.algs4.MergeX.sort(a);
        MergeX.sort(a);
        StdOut.print(timer.elapsedTime());  
        int count = 5;
        for (String x : a) {
            System.out.println(x);
            count --;
            if (count < 0) break;
        }
        
        
        Stopwatch timer2 = new Stopwatch();
        //edu.princeton.cs.algs4.Merge.sort(b);
        Merge.sort(b);
        StdOut.print(timer2.elapsedTime());  
        count = 5;
        for (String x : b) {
            System.out.println(x);
            count --;
            if (count < 0) break;
        }
        
        
    }
}
