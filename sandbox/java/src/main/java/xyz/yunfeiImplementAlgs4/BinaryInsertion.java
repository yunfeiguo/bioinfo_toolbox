
/******************************************************************************
 *  Compilation:  javac BinaryInsertion.java
 *  Execution:    java BinaryInsertion < input.txt
 *  Dependencies: StdOut.java StdIn.java
 *  Data files:   http://algs4.cs.princeton.edu/21sort/tiny.txt
 *                http://algs4.cs.princeton.edu/21sort/words3.txt
 *  
 *  Sorts a sequence of strings from standard input using 
 *  binary insertion sort with half exchanges.
 *
 *  % more tiny.txt
 *  S O R T E X A M P L E
 *
 *  % java BinaryInsertion < tiny.txt
 *  A E E L M O P R S T X                 [ one string per line ]
 *
 *  % more words3.txt
 *  bed bug dad yes zoo ... all bad yet
 *
 *  % java BinaryInsertion < words3.txt
 *  all bad bed bug dad ... yes yet zoo   [ one string per line ]
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;
import edu.princeton.cs.algs4.*;
/**
 *  The <tt>BinaryInsertion</tt> class provides a static method for sorting an
 *  array using an optimized binary insertion sort with half exchanges.
 *  <p>
 *  This implementation makes ~ N lg N compares for any array of length N.
 *  However, in the worst case, the running time is quadratic because the
 *  number of array accesses can be proportional to N^2 (e.g, if the array
 *  is reverse sorted). As such, it is not suitable for sorting large
 *  arrays (unless the number of inversions is small).
 *  <p>
 *  The sorting algorithm is stable and uses O(1) extra memory.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/21elementary">Section 2.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Ivan Pesin
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class BinaryInsertion {
    private BinaryInsertion() {}
    public static void sort(Comparable[] a) {
        if (a == null || a.length <= 1) return;
        int n = a.length;
        for (int i = 0; i < n; i++) {
            //use binary search to find the position for insertion
            //be careful about sorting stability
            //we need to consider when to stop binary search in order
            //to make the sorting stable
            //essentially, we are looking for the smallest number that
            //larger than a[i]
            int lo = 0;
            int hi = i - 1;
            while (lo <= hi) {
                int mid = lo + (hi - lo) / 2;
                if (less(a[i], a[mid])) {
                    hi = mid - 1; //a[mid] > a[i], mid could be the solution 
                } else {
                    lo = mid + 1; //a[mid] <= a[i], mid cannot be the solution
                }
            }
            for (int j = i; j > lo; j--) 
                //although number of compares is nlogn
                //number of exchanges may be quadratic in worst case (reversed sorted array)
                exch(a, j, j - 1);             
        }
    }
    /***************************************************************************
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

     // exchange a[i] and a[j]  (for indirect sort)
     private static void exch(int[] a, int i, int j) {
         int swap = a[i];
         a[i] = a[j];
         a[j] = swap;
     }
     /***************************************************************************
      *  Check if array is sorted - useful for debugging.
      ***************************************************************************/
      private static boolean isSorted(Comparable[] a) {
          return isSorted(a, 0, a.length - 1);
      }

      // is the array sorted from a[lo] to a[hi]
      private static boolean isSorted(Comparable[] a, int lo, int hi) {
          for (int i = lo+1; i <= hi; i++)
              if (less(a[i], a[i-1])) return false;
          return true;
      }

      // print array to standard output
      private static void show(Comparable[] a) {
          for (int i = 0; i < a.length; i++) {
              StdOut.println(a[i]);
          }
      }
     /**
      * Reads in a sequence of strings from standard input; insertion sorts them;
      * and prints them to standard output in ascending order.
      */
     public static void main(String[] args) {
         //String[] a = {"S", "O"};
         //String[] a = {"O", "S"};
         String[] a = {"A","O", "O","B"};
         //String[] a = {"S", "O", "R", "T", "E", "X", "A", "M", "P", "L", "E"};
         BinaryInsertion.sort(a);
         show(a);
     }

}
