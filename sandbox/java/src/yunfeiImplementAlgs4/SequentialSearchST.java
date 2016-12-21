package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac SequentialSearchST.java
 *  Execution:    java SequentialSearchST
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/31elementary/tinyST.txt  
 *  
 *  Symbol table implementation with sequential search in an
 *  unordered linked list of key-value pairs.
 *
 *  % more tinyST.txt
 *  S E A R C H E X A M P L E
 *
 *  % java SequentialSearchST < tiny.txt 
 *  L 11
 *  P 10
 *  M 9
 *  X 7
 *  H 5
 *  C 4
 *  R 3
 *  A 8
 *  E 12
 *  S 0
 *
 ******************************************************************************/

import com.sun.javafx.css.StyleCache.Key;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>SequentialSearchST</tt> class represents an (unordered)
 *  symbol table of generic key-value pairs.
 *  It supports the usual <em>put</em>, <em>get</em>, <em>contains</em>,
 *  <em>delete</em>, <em>size</em>, and <em>is-empty</em> methods.
 *  It also provides a <em>keys</em> method for iterating over all of the keys.
 *  A symbol table implements the <em>associative array</em> abstraction:
 *  when associating a value with a key that is already in the symbol table,
 *  the convention is to replace the old value with the new value.
 *  The class also uses the convention that values cannot be <tt>null</tt>. Setting the
 *  value associated with a key to <tt>null</tt> is equivalent to deleting the key
 *  from the symbol table.
 *  <p>
 *  This implementation uses a singly-linked list and sequential search.
 *  It relies on the <tt>equals()</tt> method to test whether two keys
 *  are equal. It does not call either the <tt>compareTo()</tt> or
 *  <tt>hashCode()</tt> method. 
 *  The <em>put</em> and <em>delete</em> operations take linear time; the
 *  <em>get</em> and <em>contains</em> operations takes linear time in the worst case.
 *  The <em>size</em>, and <em>is-empty</em> operations take constant time.
 *  Construction takes constant time.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/31elementary">Section 3.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
/* the use of generics is very important because it makes it possible
 * to store various types of objects as keys and values.
 */
public class SequentialSearchST<Key, Value> {
    private int N;
    private Node first;
    private class Node {
        private Key key;
        private Value val;
        private Node next;
        public Node(Key key, Value val, Node next) {
            this.key = key;
            this.val = val;
            this.next = next;
        }
    }
    public SequentialSearchST() {}
    /*
     * @return returns the number of key-value pairs in this symbol table
     */
    public int size() {
        return N;
    }
    /**    
     * @return true if this symbol table is empty
     *         false if this symbol table is not empty
     */
    public boolean isEmpty() {
        return size() == 0;
    }
    /**
     * @param key the key
     * @return the value associated with the key
     * @throws NullPointerException if key is null
     */
    public Value get(Key key) {
        if (key == null) throw new NullPointerException("argument to get() is null");
        for (Node x = first; x != null; x = x.next) {
            if (key.equals(x.key))
                return x.val;
        }
        return null;
    }
    /**
     * inserts the specified key-value pair into the symbol table, overwriting any existing value
     * associated the key, deletes the key if the value==null
     * @param key they key
     * @param value the value
     * @throws NullPointerException
     */
    public void put(Key key, Value value) {
        if (key == null) throw new NullPointerException("first argument must not be null");
        if (value == null) {
            delete(key);
            return;
        }
        for (Node x = first; x != null; x = x.next) {
            if (key.equals(x.key)) {
                x.val = value;
                return;
            }
        }
        first = new Node(key, value, first);
        N++;
    }
    /**
     * removes the specified key and its associated value from this symbol table
     * if they exist
     * @param key
     * @throws NullPointerException if <tt>key</tt> is <tt>null</tt>
     */
    public void delete(Key key) {
        if (key == null) throw new NullPointerException();
        first = delete(first, key);        
    }
    /**
     * removes the specified key from the list with first as its head
     * @param first head of list
     * @param key the key to be deleted
     * @return new head
     */
    private Node delete(Node first, Key key) {
        if (key == null) return null;
        Node anchor = new Node(null, null, first);
        Node pointer = anchor;
        while (pointer.next != null) {
            if (key.equals(pointer.next.key)) {
                pointer.next = pointer.next.next;
                N--;
                break;
            }
            pointer = pointer.next;
        }
        return anchor.next;
    }
    /**
     * returns all keys in the symbol table as an <tt>Iterable</tt>.
     * 
     * @return all keys in the symbol table
     */
    public Iterable<Key> keys() {
        Queue<Key> queue = new Queue<Key>();
        for (Node x = first; x != null; x = x.next)
            queue.enqueue(x.key);
        return queue;
    }
    public static void main(String[] args) {
        SequentialSearchST<String, Integer> st = new SequentialSearchST<>();
        for (int i = 0; !StdIn.isEmpty(); i++) {
            String key = StdIn.readString();
            st.put(key, i);
        }
        for (String s: st.keys()) {
            StdOut.println(s + " " + st.get(s));
        }
    }
}
