/******************************************************************************
 *  Compilation:  javac Shell.java
 *  Execution:    java Shell < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/21sort/tiny.txt
 *                http://algs4.cs.princeton.edu/21sort/words3.txt
 *   
 *  Sorts a sequence of strings from standard input using shellsort.
 *
 *  Uses increment sequence proposed by Sedgewick and Incerpi.
 *  The nth element of the sequence is the smallest integer >= 2.5^n
 *  that is relatively prime to all previous terms in the sequence.
 *  For example, incs[4] is 41 because 2.5^4 = 39.0625 and 41 is
 *  the next integer that is relatively prime to 3, 7, and 16.
 *   
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java Shell < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *    
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *  
 *  % java Shell < words3.txt
 *  all bad bed bug dad ... yes yet zoo    [ one string per line ]
 *
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdIn;

/**
 *  The <tt>Shell</tt> class provides static methods for sorting an
 *  array using Shellsort with Knuth's increment sequence (1, 4, 13, 40, ...).
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Shell {
    public static void main(String[] args) {
        String[] a = StdIn.readAllStrings();
        Shell.sort(a);
        show(a);
    }
    private static void show(Comparable[] a) {
        for (Comparable i : a) 
            System.out.println(i);
    }
    
    private static boolean isHsorted(Comparable[] a, int h) {
        for (int i = h; i < a.length; i++)
            if (less(a[i], a[i - h]))
                return false;
        return true;
    }
    @SuppressWarnings("unchecked")
    public static boolean less(Comparable a, Comparable b) {
        return a.compareTo(b) < 0;
    }
    
    public static void exch(Comparable[] a, int x, int y) {
        Comparable tmp = a[x];
        a[x] = a[y];
        a[y] = tmp;
    }
    public static void sort(Comparable[] a) {
        int h = 1;
        while (h < a.length/3) h = h * 3 + 1; //get h
        while (h >= 1) {
            for (int i = h; i < a.length; i++) {
                for (int j = i; j >= h && less(a[j], a[j - h]); j -= h)
                    exch(a, j, j - h);
            }
            h /= 3;
        }            
    }

}
