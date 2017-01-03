package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac LinkedBag.java
 *  Execution:    java LinkedBag < input.txt
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


import java.security.NoSuchAlgorithmException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
/**
 *  The <tt>LinkedBag</tt> class represents a bag (or multiset) of 
 *  generic items. It supports insertion and iterating over the 
 *  items in arbitrary order.
 *  <p>
 *  This implementation uses a singly-linked list with a non-static nested class Node.
 *  See {@link Bag} for a version that uses a static nested class.
 *  The <em>add</em>, <em>isEmpty</em>, and <em>size</em> operations
 *  take constant time. Iteration takes time proportional to the number of items.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/13stacks">Section 1.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class LinkedBag<Item> implements Iterable<Item> {
        private int N;
        private Node head;
        
        private class Node {
            private Item val;
            private Node next;
        }
        
        public LinkedBag() {
            N = 0;
            head = null;
        }
        
        //add Item to head of linked list
        public void add(Item x) {
            Node newOne = new Node();
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
        }
        
        private class ListIterator implements Iterator<Item> {
            private Node cur  = head;
            
            public ListIterator() {
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
        public static void main(String[] args) {
            LinkedBag<String> bag = new LinkedBag<String>();
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
