/******************************************************************************
 *  Compilation:  javac Selection.java
 *  Execution:    java  Selection < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/21sort/tiny.txt
 *                http://algs4.cs.princeton.edu/21sort/words3.txt
 *   
 *  Sorts a sequence of strings from standard input using selection sort.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java Selection < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *    
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *  
 *  % java Selection < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;

import java.util.Comparator;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

public class Selection {
    /**
     *  The <tt>Selection</tt> class provides static methods for sorting an
     *  array using selection sort.
     *  <p>
     *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *
     *  @author Robert Sedgewick
     *  @author Kevin Wayne
     */
    //public Selection() {}
    private Selection() {} //all methods are static, this class should not be instantiated
    @SuppressWarnings("unchecked")
    public static boolean less(Comparable a, Comparable b) {
        return a.compareTo(b) < 0;
    }
    public static boolean less(Comparator c, Object a, Object b) {
        return c.compare(a, b) < 0;
    }
    // print array to standard output
    private static void show(Comparable[] a) {
        for (int i = 0; i < a.length; i++) {
            StdOut.println(a[i]);
        }
    }
    public static void exch(Object[] a, int x, int y) {
        Object tmp = a[x];
        a[x] = a[y];
        a[y] = tmp;
    }
    public static void sort(Comparable[] a) {
        for (int i = 0; i < a.length; i++) {
            //each time, we find the min to the right (inclusive) of i
            //put it to it, then move forward
            int min = i;
            for (int j = i + 1; j < a.length; j++)
                if (less(a[j], a[min])) min = j;
            exch(a, i, min);
        }
    }
    public static void sort(Object[] a, Comparator c) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            int min = i;
            for (int j = i + 1; j < N; j++)
                if (less(c, a[j], a[min])) min = j;
            exch(a, i, min);
        }
    }
    public static void main(String[] args) {
        Selection test = new Selection();
        String[] a = StdIn.readAllStrings();
        Selection.sort(a);
        show(a);

    }
}
