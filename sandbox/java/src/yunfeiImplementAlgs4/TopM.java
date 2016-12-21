package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac TopM.java
 *  Execution:    java TopM M < input.txt
 *  Dependencies: MinPQ.java Transaction.java StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/24pq/tinyBatch.txt
 * 
 *  Given an integer M from the command line and an input stream where
 *  each line contains a String and a long value, this MinPQ client
 *  prints the M lines whose numbers are the highest.
 * 
 *  % java TopM 5 < tinyBatch.txt 
 *  Thompson    2/27/2000  4747.08
 *  vonNeumann  2/12/1994  4732.35
 *  vonNeumann  1/11/1999  4409.74
 *  Hoare       8/18/1992  4381.21
 *  vonNeumann  3/26/2002  4121.85
 *
 ******************************************************************************/

import java.util.Deque;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Stack;

import edu.princeton.cs.algs4.MinPQ;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Transaction;

/**
 *  The <tt>TopM</tt> class provides a client that reads a sequence of
 *  transactions from standard input and prints the <em>M</em> largest ones
 *  to standard output. This implementation uses a {@link MinPQ} of size
 *  at most <em>M</em> + 1 to identify the <em>M</em> largest transactions
 *  and a {@link Stack} to output them in the proper order.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/24pq">Section 2.4</a>
 *  of <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class TopM {
    private TopM() {}
    /**
     *  Reads a sequence of transactions from standard input; takes a
     *  command-line integer M; prints to standard output the M largest
     *  transactions in descending order.
     */
    /*public static void main(String[] args) {
        int M = Integer.parseInt(args[0]);
        MinPQ<Transaction> pq = new MinPQ<Transaction>(M+1);
        while (StdIn.hasNextLine()) {
            String line = StdIn.readLine();
            Transaction transaction = new Transaction(line);
            pq.insert(transaction);
            if (pq.size() > M) {
                pq.delMin();
            }
        }
        Stack<Transaction> stack = new Stack<Transaction>();
        for (Transaction transaction : pq) {
            stack.push(transaction);            
        }
        for (Transaction transaction : stack)
            StdOut.println(transaction);
    }*/
    //implementation using Java builtin priority queue
    public static void main(String[] args) {
        int M = Integer.parseInt(args[0]);
        PriorityQueue<Transaction> pq = new PriorityQueue<Transaction>(M+1);
        while (StdIn.hasNextLine()) {
            String line = StdIn.readLine();
            Transaction transaction = new Transaction(line);
            pq.offer(transaction);
            if (pq.size() > M) {
                pq.poll();
            }
        }
        Deque<Transaction> stack = new LinkedList<Transaction>();
        //Java implementation of priority queue iterator does not have order
        for (Transaction transaction : pq) {
            stack.addFirst(transaction);            
        }
        for (Transaction transaction : stack)
            StdOut.println(transaction);
    }

}
