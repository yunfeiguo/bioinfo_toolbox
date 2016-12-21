package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

public class Insertion {
    private Insertion() { }
    public static boolean less(Comparable a, Comparable b) {
        return a.compareTo(b) < 0;
    }
    
    public static void exch(Object[] a, int x, int y) {
        Object tmp = a[x];
        a[x] = a[y];
        a[y] = tmp;
    }
    public static void sort(Comparable[] a) {
        int N = a.length;
        for (int i = 0; i < N; i++) {
            //each iteration we make sure
            //elements 0~i are in order
            //by moving small elements to the leftmost
            //position allowed
            for (int j = i; j > 0; j--)
                if (less(a[j], a[j - 1]))
                    exch(a, j, j - 1);  
                else
                    break;
        }
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
        String[] a = StdIn.readAllStrings();
        Insertion.sort(a);
        //Insertion.sort(null); //not handled
        show(a);
    }
}
