package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac BinarySearchST.java
 *  Execution:    java BinarySearchST
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/31elementary/tinyST.txt  
 *  
 *  Symbol table implementation with binary search in an ordered array.
 *
 *  % more tinyST.txt
 *  S E A R C H E X A M P L E
 *  
 *  % java BinarySearchST < tinyST.txt
 *  A 8
 *  C 4
 *  E 12
 *  H 5
 *  L 11
 *  M 9
 *  P 10
 *  R 3
 *  S 0
 *  X 7
 *
 ******************************************************************************/

import java.util.NoSuchElementException;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>BST</tt> class represents an ordered symbol table of generic
 *  key-value pairs.
 *  It supports the usual <em>put</em>, <em>get</em>, <em>contains</em>,
 *  <em>delete</em>, <em>size</em>, and <em>is-empty</em> methods.
 *  It also provides ordered methods for finding the <em>minimum</em>,
 *  <em>maximum</em>, <em>floor</em>, <em>select</em>, and <em>ceiling</em>.
 *  It also provides a <em>keys</em> method for iterating over all of the keys.
 *  A symbol table implements the <em>associative array</em> abstraction:
 *  when associating a value with a key that is already in the symbol table,
 *  the convention is to replace the old value with the new value.
 *  Unlike {@link java.util.Map}, this class uses the convention that
 *  values cannot be <tt>null</tt>&mdash;setting the
 *  value associated with a key to <tt>null</tt> is equivalent to deleting the key
 *  from the symbol table.
 *  <p>
 *  This implementation uses a sorted array. It requires that
 *  the key type implements the <tt>Comparable</tt> interface and calls the
 *  <tt>compareTo()</tt> and method to compare two keys. It does not call either
 *  <tt>equals()</tt> or <tt>hashCode()</tt>.
 *  The <em>put/em> and <em>remove</em> operations each take linear time in
 *  the worst case; the <em>contains</em>, <em>ceiling</em>, <em>floor</em>,
 *  and <em>rank</em> operations take logarithmic time; the <em>size</em>,
 *  <em>is-empty</em>, <em>minimum</em>, <em>maximum</em>, and <em>select</em>
 *  operations take constant time. Construction takes constant time.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/31elementary">Section 3.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  For other implementations, see {@link ST}, {@link BST},
 *  {@link SequentialSearchST}, {@link RedBlackBST},
 *  {@link SeparateChainingHashST}, and {@link LinearProbingHashST},
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */
public class BinarySearchST<Key extends Comparable<Key>, Value> {
    private static final int INIT_CAPACITY = 2;
    private Key[] keys;
    private Value[] vals;
    private int N = 0;
    public BinarySearchST() {
        this(INIT_CAPACITY);
    }
    public BinarySearchST(int capacity) {
        keys = (Key[]) new Comparable[capacity];
        vals = (Value[]) new Object[capacity];
    }
    private void resize(int capacity) {
        assert capacity >= N;
        Key[] newKey = (Key[]) new Comparable[capacity];
        Value[] newVal = (Value[]) new Object[capacity];
        for (int i = 0; i < N; i++) {
            newKey[i] = keys[i];
            newVal[i] = vals[i];
        }
        keys = newKey;
        vals = newVal;
        return;
    }
    public int size() {
        return N;
    }
    public boolean isEmpty() {
        return size() == 0;
    }
    public boolean contains(Key key) {
        if (key == null) throw new NullPointerException("argument to contains() is null");
        return get(key) != null;
    }
    public Value get(Key key) {
        if (key == null) throw new NullPointerException("argument cannot be null");
        if (isEmpty()) return null;
        int i = rank(key);
        if (i < N && keys[i].compareTo(key) == 0)
            return vals[i];
        return null;
    }
    public int rank(Key key) {
        if (key == null) throw new NullPointerException();
        int lo = 0, hi = N - 1;
        while (lo <= hi) {
            int mid = lo + (hi - lo)/2;
            int comp = key.compareTo(keys[mid]);
            if (comp < 0) hi = mid - 1;
            else if (comp > 0) lo = mid + 1;
            else return mid;
        }
        return lo; //we want to make sure return an index
                   //where all elements smaller to key appear
                   //in the left part
    }
    public void put(Key key, Value val) {
        if (key == null) throw new NullPointerException("put() cannot be null");
        if (val == null) {
            delete(key);
            return;
        }
        int i = rank(key);
        if (i < N && keys[i].compareTo(key) == 0) {
            vals[i] = val;
            return;
        }
        if (N == keys.length) {
            resize(2*N);
        }
        for (int j = N; j > i; j--) {
            keys[j] = keys[j - 1];
            vals[j] = vals[j - 1];
        }
        keys[i] = key;
        vals[i] = val;
        N++;
        assert(check());
        return;
    }
    public void delete(Key key) {
        if (key == null) throw new NullPointerException("delete() null");
        if (isEmpty()) return;
        int i = rank(key);        
        if (i < N && keys[i].compareTo(key) == 0) {            
            for (int j = i; j < N - 1; j++) {
                keys[j] = keys[j+1];
                vals[j] = vals[j+1];
            }
            keys[N] = null;
            vals[N] = null;
            N--;
            if (N > 0 && N <= keys.length/4) {
                resize(keys.length/2);
            }
        } else {
            //didn't find the key, do nothing
            return;
        }
    }
    public void deleteMin() {
        if (isEmpty()) throw new NoSuchElementException();
        delete(min());
    }
    public void deleteMax() {
        if (isEmpty()) throw new NoSuchElementException();
        delete(max());
    }
    public Key min() {
        if (isEmpty()) return null;
        return keys[0];
    }
    public Key max() {
        if (isEmpty()) return null;
        return keys[N-1];
    }
    public Key select(int k) {
        if (k < 0 || k >= N) return null;
        return keys[k];
    }
    public Key floor(Key key) {
        if (key == null) throw new NullPointerException();
        int i = rank(key);
        if (i < N && keys[i].compareTo(key) == 0) return keys[i];
        if (i == 0) return null;
        else return keys[i - 1]; //notice here we always return the largest element smaller than key        
    }
    public Key ceiling(Key key) {
        if (key == null) throw new NullPointerException();
        int i = rank(key);
        //think about it: why don't we return keys[i+1]???
        if (i == N) return null;
        else return keys[i];
    }
    public int size(Key lo, Key hi) {
        if (lo == null || hi == null) throw new NullPointerException();
        if (lo.compareTo(hi) > 0) return 0;
        if (contains(hi)) return rank(hi) - rank(lo) + 1;
        else return rank(hi) - rank(lo);
    }
    public Iterable<Key> keys() {
        return keys(min(), max());
    }
    public Iterable<Key> keys(Key lo, Key hi) {
        if (lo == null || hi == null) throw new NullPointerException();
        Queue<Key> queue = new Queue<Key>();
        if (lo.compareTo(hi) > 0) return queue;        
        for (int i = rank(lo); i < N && keys[i].compareTo(hi) <= 0; i++) {
            queue.enqueue(keys[i]);
        }
        return queue;
    }
    /***************************************************************************
     *  Check internal invariants.
     ***************************************************************************/

     private boolean check() {
         return isSorted() && rankCheck();
     }

     // are the items in the array in ascending order?
     private boolean isSorted() {
         for (int i = 1; i < size(); i++)
             if (keys[i].compareTo(keys[i-1]) < 0) return false;
         return true;
     }

     // check that rank(select(i)) = i
     private boolean rankCheck() {
         for (int i = 0; i < size(); i++)
             if (i != rank(select(i))) return false;
         for (int i = 0; i < size(); i++)
             if (keys[i].compareTo(select(rank(keys[i]))) != 0) return false;
         return true;
     }
     /**
      * Unit tests the <tt>BinarySearchST</tt> data type.
      */
     public static void main(String[] args) { 
         BinarySearchST<String, Integer> st = new BinarySearchST<String, Integer>();
         for (int i = 0; !StdIn.isEmpty(); i++) {
             String key = StdIn.readString();
             st.put(key, i);
         }
         for (String s : st.keys())
             StdOut.println(s + " " + st.get(s));
     }

}
