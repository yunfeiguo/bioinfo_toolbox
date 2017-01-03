package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac Heap.java
 *  Execution:    java Heap < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/24pq/tiny.txt
 *                http://algs4.cs.princeton.edu/24pq/words3.txt
 *  
 *  Sorts a sequence of strings from standard input using heapsort.
 *
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java Heap < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *
 *  % more words3.txt
 * bed bug dad yes zoo
 * now for tip ilk dim 
 *tag jot sob nob sky
 *hut men egg few jay
 *owl joy rap gig wee
 *was wad fee tap tar
 *dug jam all bad yet
 *
 *  % java Heap < words3.txt
 *  all bad bed bug dad ... yes yet zoo   [ one string per line ]
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>Heap</tt> class provides a static methods for heapsorting
 *  an array.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/24pq">Section 2.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Heap {
    // This class should not be instantiated.
    private Heap() { }

    /**
     * Rearranges the array in ascending order, using the natural order.
     * @param pq the array to be sorted
     */
    public static void sort(Comparable[] pq) {
        //first establish heap order
        //then take the max, put it at the end
        int parent = pq.length/2;
        //heapify the array
        while (parent > 0) {
            sink(pq, parent, pq.length);
            parent--;
        }
        for (int i = pq.length;i > 1; i--) {
            exch(pq, 1, i);
            sink(pq, 1, i - 1);
        }
        assert(isSorted(pq));
    }
    /***************************************************************************
     * Helper functions to restore the heap invariant.
     ***************************************************************************/
    private static void sink(Comparable[] pq, int k, int N) {
        //bottom up sink to restore heap order, k is root, N is length
        //max heap
        while (2*k <= N) {
            int child = 2*k;
            if (child < N && less(pq, child, child + 1))
                child++;
            if (!less(pq, k, child)) //there are two stopping criteria, one is about index, the other about order
                break;
            exch(pq, child, k);
            k = child;
        }
    }
    /***************************************************************************
     * Helper functions for comparisons and swaps.
     * Indices are "off-by-one" to support 1-based indexing.
     ***************************************************************************/
    private static boolean less(Comparable[] pq, int i, int j) {
        return less(pq[i - 1], pq[j - 1]);
    }
    private static void exch(Object[] pq, int i, int j) {
        Object tmp = pq[i - 1];
        pq[i - 1] = pq[j - 1];
        pq[j- 1] = tmp;
    }
    // is v < w ?
    private static boolean less(Comparable v, Comparable w) {
               return v.compareTo(w) < 0;
    }
    /***************************************************************************
     *  Check if array is sorted - useful for debugging.
     ***************************************************************************/
    private static boolean isSorted(Comparable[] a) {
        if (a == null || a.length == 1) {
            return true;
        }
        for (int i = 1; i < a.length; i++) {
            if (less(a[i],a[i - 1])) {
                return false;
            }            
        }
        return true;
    }

    // print array to standard output
    private static void show(Comparable[] a) {
        for (Comparable i : a) {
            System.out.println(i);
        }
    }
    /**
     * Reads in a sequence of strings from standard input; heapsorts them; 
     * and prints them to standard output in ascending order. 
     */
    public static void main(String[] args) {
        String[] a = StdIn.readAllStrings();
        Heap.sort(a);
        show(a);
    }
}
