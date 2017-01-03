package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac MinPQ.java
 *  Execution:    java MinPQ < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  
 *  Generic min priority queue implementation with a binary heap.
 *  Can be used with a comparator instead of the natural order.
 *  % more tinyPQ.txt
 *  P Q E - X A M - P L E -
 *  
 *  % java MinPQ < tinyPQ.txt
 *  E A E (6 left on pq)
 *
 *  We use a one-based array to simplify parent and child calculations.
 *
 *  Can be optimized by replacing full exchanges with half exchanges
 *  (ala insertion sort).
 *
 ******************************************************************************/

import java.util.Comparator;
import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 *  The <tt>MinPQ</tt> class represents a priority queue of generic keys.
 *  It supports the usual <em>insert</em> and <em>delete-the-minimum</em>
 *  operations, along with methods for peeking at the minimum key,
 *  testing if the priority queue is empty, and iterating through
 *  the keys.
 *  <p>
 *  This implementation uses a binary heap.
 *  The <em>insert</em> and <em>delete-the-minimum</em> operations take
 *  logarithmic amortized time.
 *  The <em>min</em>, <em>size</em>, and <em>is-empty</em> operations take constant time.
 *  Construction takes time proportional to the specified capacity or the number of
 *  items used to initialize the data structure.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/24pq">Section 2.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 *
 *  @param <Key> the generic type of key on this priority queue
 */
public class MinPQ<Key extends Comparable<Key>> implements Iterable<Key> {
    private Key[] pq; //one-based array for priority queue
    private int N;
    private Comparator<Key> comparator;
    /**
     * Initializes an empty priority queue with the given initial capacity.
     *
     * @param  initCapacity the initial capacity of this priority queue
     */
    @SuppressWarnings("unchecked")
    public MinPQ(int initCapacity) {
        pq = (Key[]) new Comparable[initCapacity + 1];
        N = 0;
    }
    /**
     * Initializes an empty priority queue.
     */
    @SuppressWarnings("unchecked")
    public MinPQ() {
        pq = (Key[]) new Comparable[2];
        N = 0;
    }       
    /**
     * Initializes an empty priority queue with the given initial capacity,
     * using the given comparator.
     *
     * @param  initCapacity the initial capacity of this priority queue
     * @param  comparator the order to use when comparing keys
     */
    @SuppressWarnings("unchecked")
    public MinPQ(int initCapacity, Comparator<Key> comparator) {
        pq = (Key[]) new Comparable[initCapacity + 1];
        N = 0;
        this.comparator = comparator;
    }        

    /**
     * Initializes an empty priority queue using the given comparator.
     *
     * @param  comparator the order to use when comparing keys
     */
    @SuppressWarnings("unchecked")
    public MinPQ(Comparator<Key> comparator) {
        pq = (Key[]) new Comparable[2];
        N = 0;
        this.comparator = comparator;
    }

    /**
     * Initializes a priority queue from the array of keys.
     * <p>
     * Takes time proportional to the number of keys, using sink-based heap construction.
     *
     * @param  keys the array of keys
     */
    @SuppressWarnings("unchecked")
    public MinPQ(Key[] keys) {
        pq = (Key[]) new Comparable[keys.length + 1];
        /* this is one way for construction
         * however it takes O(NlgN) time
         */
        /*
        for (Key k : keys) {
            pq.insert(k);
        }
        */
        for (int i = 0; i < keys.length; i++) {
            pq[i] = keys[i];
        }
        N = keys.length;
        //we begin with last node with a child
        //make the tree in heap order
        //then gradually go up
        //traverse each level bottom up
        for (int i = N/2; i >= 1; i--) {
            sink(i);
        }           
        assert isMinHeap();
        //time complexity?
        //suppose N = 2^k
        //# of exchanges = 1*2^(k-2) + 2*2^(k-3) + ... + k*2^(0)
        //multiply it by 2, do subtraction of original
        //# of exch = N - k < N
        //# of compares is 2 times as large as # of exch.
        //space complexity O(N)
    }
    /**
     * Returns true if this priority queue is empty.
     *
     * @return <tt>true</tt> if this priority queue is empty;
     *         <tt>false</tt> otherwise
     */
    public boolean isEmpty() {
        return N == 0;
    }

    /**
     * Returns the number of keys on this priority queue.
     *
     * @return the number of keys on this priority queue
     */
    public int size() {
        return N;
    }
    /**
     * Returns a smallest key on this priority queue.
     *
     * @return a smallest key on this priority queue
     * @throws NoSuchElementException if this priority queue is empty
     */
    public Key min() {
        if (isEmpty())
            throw new NoSuchElementException();
        return pq[1];
    }
    // helper function to double the size of the heap array
    private void resize(int capacity) {
        Key[] pq2 = (Key[]) new Comparable[capacity + 1];
        for (int i = 1; i < pq.length; i++) {
            pq2[i] = pq[i];
        }
        pq = pq2;        
    }
    /**
     * Adds a new key to this priority queue.
     *
     * @param  x the key to add to this priority queue
     */
    public void insert(Key x) {
        if (N == pq.length - 1) resize(N*2);
        pq[++N] = x;
        swim(N);
        assert isMinHeap();
    }
    /**
     * Removes and returns a smallest key on this priority queue.
     *
     * @return a smallest key on this priority queue
     * @throws NoSuchElementException if this priority queue is empty
     */
    public Key delMin() {
        Key result = pq[1];
        exch(1, N);
        pq[N--] = null; //avoid loitering        
        sink(1);
        if (N > 0 && N < (pq.length - 1)/4) resize(N/2);
        assert isMinHeap();
        return result;
    }
    /***************************************************************************
     * Helper functions to restore the heap invariant.
     ***************************************************************************/

     private void swim(int k) {
         //swim up
         assert(k <= N && k >= 1);
         while (k <= N && k >= 1 && k/2 >= 1 && greater(k/2, k)) {
             exch(k/2, k);
             k = k/2;
         }
     }
     private void sink(int k) {
         //sink down
         assert(k >= 1 && k <= N);
         while (k >= 1 && k <= N && k*2 <= N) {
             int childIndex = k*2;
             if (childIndex < N && greater(childIndex, childIndex + 1)) {
                 //exchange with smaller child
                 childIndex++;
             }
             if (greater(k, childIndex)) {
                 exch(k, childIndex);
                 k = childIndex;
             } else {
                 break;
             }
         }
     }
     /***************************************************************************
      * Helper functions for compares and swaps.
      ***************************************************************************/
      private boolean greater(int i, int j) {
          if (comparator != null) {
              return comparator.compare(pq[i], pq[j]) > 0;
          } else {
              return pq[i].compareTo(pq[j]) > 0;
          }
      }
      private void exch(int i, int j) {
          Key tmp = pq[i];
          pq[i] = pq[j];
          pq[j] = tmp;
      }
      // is pq[1..N] a min heap?
      private boolean isMinHeap() {
          System.out.println("isMinHeap(1) " + isMinHeap(1));
          return isMinHeap(1);
      }

      // is subtree of pq[1..N] rooted at k a min heap?
      private boolean isMinHeap(int k) {
          /*
          assert(k <= N && k >= 1);
          if (k <= N && 2*k > N) {
              //k is leaf node
              return true;
          } else {
              //only left child
              if (2*k == N && !greater(k, k*2)) {
                  return true;
              }
              //both children are present
              if (2*k < N && !greater(k, k*2) && !greater(k, k*2 + 1)) {
                  return isMinHeap(2*k) && isMinHeap(2*k + 1);
              }
              return false;
          }*/
          if (k > N) return true;
          int left = 2*k;
          int right = 2*k + 1;
          if (left <= N && greater(k, left)) return false;
          if (right <= N && greater(k, right)) return false;
          return isMinHeap(left) && isMinHeap(right);
      }
      /**
       * Returns an iterator that iterates over the keys on this priority queue
       * in ascending order.
       * <p>
       * The iterator doesn't implement <tt>remove()</tt> since it's optional.
       *
       * @return an iterator that iterates over the keys in ascending order
       */
      public Iterator<Key> iterator() {
          return new HeapIterator();
      }
      private class HeapIterator implements Iterator<Key> {
          private MinPQ<Key> iteratorPQ;
          public HeapIterator() {
              iteratorPQ = new MinPQ<Key>(N);
              for (int i = 1; i <= N; i++) {
                  iteratorPQ.insert(pq[i]);
              }
          }
          public void remove() {throw new NotImplementedException();}
          public boolean hasNext() { return !iteratorPQ.isEmpty(); }
          @SuppressWarnings("unchecked")
        public Key next() {
              if(!hasNext()) {throw new NoSuchElementException();}
              return (Key) iteratorPQ.delMin();
          }
      }
      /**
       * Unit tests the <tt>MinPQ</tt> data type.
       */
      public static void main(String[] args) {
          MinPQ<String> pq = new MinPQ<String>();
          while (!StdIn.isEmpty()) {
              String item = StdIn.readString();
              if (!item.equals("-")) pq.insert(item);
              else if (!pq.isEmpty()) StdOut.print(pq.delMin() + " ");
          }
          StdOut.println("(" + pq.size() + " left on pq)");
      }
}
