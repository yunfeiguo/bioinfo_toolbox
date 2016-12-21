package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac UnorderedArrayMaxPQ.java
 *  Execution:    java UnorderedArrayMaxPQ
 *  Dependencies: StdOut.java 
 *  
 *  Priority queue implementation with an unsorted array.
 * 
 *  Limitations
 *  -----------
 *   - no array resizing
 *   - does not check for overflow or underflow.
 *
 ******************************************************************************/
public class UnorderedArrayMaxPQ<Key extends Comparable<Key>> {
    private Key[] pq;
    private int N;
    
    @SuppressWarnings("unchecked")
    public UnorderedArrayMaxPQ(int capacity) {
        pq = (Key[]) new Comparable[capacity];
        N = 0;
    }
    public boolean isEmpty() { return N == 0; }
    public int size() { return N; }
    public void insert(Key x) { pq[N++] = x; } //constant time
    //removing all elements takes ~N^2/2 operations
    public Key delMax() {
        //this operation takes ~N operations
        int max = 0;
        if (isEmpty()) ;//throw exception
        for (int i = 1; i < N; i++) {
            if (less(max, i)) {
                max = i;
            }
        }
        exch(max, N - 1);
        Key result = pq[N - 1];
        pq[N - 1] = null; //avoid loitering
        N--;
        return result;              
    }
    /***************************************************************************
     * Helper functions.
     ***************************************************************************/
     private boolean less(int i, int j) {
         return pq[i].compareTo(pq[j]) < 0;
     }

     private void exch(int i, int j) {
         Key swap = pq[i];
         pq[i] = pq[j];
         pq[j] = swap;
     }
     /***************************************************************************
      * Test routine.
      ***************************************************************************/
      public static void main(String[] args) {
          UnorderedArrayMaxPQ<String> pq = new UnorderedArrayMaxPQ<String>(10);
          pq.insert("this");
          pq.insert("is");
          pq.insert("a");
          pq.insert("test");
          while (!pq.isEmpty()) 
              StdOut.println(pq.delMax());
      }
}
