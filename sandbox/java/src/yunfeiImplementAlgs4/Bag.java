/******************************************************************************
 *  Compilation:  javac Bag.java
 *  Execution:    java Bag < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *
 *  A generic bag or multiset, implemented using a singly-linked list.
 *
 *  % more tobe.txt 
 *  to be or not to - be - - that - - - is
 *
 *  % java Bag < tobe.txt
 *  size of bag = 14
 *  is
 *  -
 *  -
 *  -
 *  that
 *  -
 *  -
 *  be
 *  -
 *  to
 *  not
 *  or
 *  be
 *  to
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;

import java.security.NoSuchAlgorithmException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>Bag</tt> class represents a bag (or multiset) of 
 *  generic items. It supports insertion and iterating over the 
 *  items in arbitrary order.
 *  <p>
 *  This implementation uses a singly-linked list with a static nested class Node.
 *  See {@link LinkedBag} for the version from the
 *  textbook that uses a non-static nested class.
 *  The <em>add</em>, <em>isEmpty</em>, and <em>size</em> operations
 *  take constant time. Iteration takes time proportional to the number of items.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/13stacks">Section 1.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 *
 *  @param <Item> the generic type of an item in this bag
 */
public class Bag<Item> implements Iterable<Item> {
    private int N;
    private Node<Item> head;
    
    private static class Node<Item> {
        private Item val;
        private Node<Item> next;
    }
    
    public Bag() {
        N = 0;
        head = null;
    }
    
    //add Item to head of linked list
    public void add(Item x) {
        Node<Item> newOne = new Node<Item>();
        newOne.val = x;
        newOne.next = head;
        head = newOne;
        N++;
    }
    
    public boolean isEmpty() {
        return N == 0;
    }
    
    public int size() {
        return N;
    }
    
    public Iterator<Item> iterator() {
        return new ListIterator();
        //return new StaticListIterator(head);
    }
    /*
     * static nested class has no access to anything in the enclosing class, e.g. the Item
    private static class StaticListIterator implements Iterator<Item> {
        private Node<Item> cur;
        
        public StaticListIterator(Node<Item> head) {
            cur = head;
        }
        public boolean hasNext() {
            return cur != null;
        }
        public Item next() {
            if (!hasNext()) throw new NoSuchElementException();
            Item ret = cur.val;
            cur = cur.next;
            return ret;
        }
        public void remove() {
            throw new UnsupportedOperationException();
        }
        
    }
    */
    
    private class ListIterator implements Iterator<Item> {
        private Node<Item> cur = head;

        public boolean hasNext() {
            return cur != null;
        }
        public Item next() {
            if (!hasNext()) throw new NoSuchElementException();
            Item ret = cur.val;
            cur = cur.next;
            return ret;
        }
        public void remove() {
            throw new UnsupportedOperationException();
        }
        
    }
    
    public static void main(String[] args) {
        Bag<String> bag = new Bag<String>();
        while (!StdIn.isEmpty()) {
            String item = StdIn.readString();
            bag.add(item);
        }
        StdOut.println("size of bag = " + bag.size());
        for (String s : bag) {
            StdOut.println(s);
        }
    }

}
