/******************************************************************************
 *  Compilation:  javac Queue.java
 *  Execution:    java Queue < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/13stacks/tobe.txt  
 *
 *  A generic queue, implemented using a linked list.
 *
 *  % java Queue < tobe.txt 
 *  to be or not to be (2 left on queue)
 *
 ******************************************************************************/

package yunfeiImplementAlgs4;

import java.util.Iterator;
import java.util.NoSuchElementException;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>Queue</tt> class represents a first-in-first-out (FIFO)
 *  queue of generic items.
 *  It supports the usual <em>enqueue</em> and <em>dequeue</em>
 *  operations, along with methods for peeking at the first item,
 *  testing if the queue is empty, and iterating through
 *  the items in FIFO order.
 *  <p>
 *  This implementation uses a singly-linked list with a static nested class for
 *  linked-list nodes. See {@link LinkedQueue} for the version from the
 *  textbook that uses a non-static nested class.
 *  The <em>enqueue</em>, <em>dequeue</em>, <em>peek</em>, <em>size</em>, and <em>is-empty</em>
 *  operations all take constant time in the worst case.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/13stacks">Section 1.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 *
 *  @param <Item> the generic type of an item in this bag
 */
public class Queue<Item> implements Iterable<Item> {
    private int N;
    private Node<Item> head;
    private Node<Item> tail;
    
    private static class Node<Item> { //here we must specify <Item> generic, otherwise we can't use Item as a field
        private Item val;
        private Node<Item> next;
    }
    
    public Queue() {
        N = 0;
        head = null;
        tail = null;
    }
    
    public boolean isEmpty() {
        return N == 0;
    }
    
    public int size() {
        return N;
    }
    
    public Item peek() {
        if (isEmpty()) throw new NoSuchElementException();
        return tail.val;
    }
    
    public void enqueue(Item item) {
        Node<Item> newNode = new Node<Item>();
        newNode.val = item;
        newNode.next = null;
        if (isEmpty()) head = tail = newNode;
        else {
            tail.next = newNode;
            tail = newNode;
        }
        N++;
    }
    
    public Item dequeue() {
        if (isEmpty()) throw new NoSuchElementException();
        Item ret = head.val;
        head = head.next;
        N--;
        return ret;
    }
    
    public Iterator<Item> iterator() { //only class can 'IMPLEMENT' interface
        return new ListIterator<Item>(head);
        //return new ListIterator<Item>();
    }
    
    private class ListIterator<Item> implements Iterator<Item> {
        private Node<Item> cur;
        //private Node<Item> cur = head; //this also works, but with less encapsulation
        
        public ListIterator(Node<Item> head) {
            cur = head;
        }     
        
        public boolean hasNext() { return cur != null; }
        public void remove() { throw new UnsupportedOperationException(); }
        public Item next() {
            if (!hasNext()) throw new NoSuchElementException();
            Item ret = cur.val;
            cur = cur.next;
            return ret;
        }
    }
    /**
     * Unit tests the <tt>Queue</tt> data type.
     */
    public static void main(String[] args) {
        Queue<String> q = new Queue<String>();
        while (!StdIn.isEmpty()) {
            String item = StdIn.readString();
            if (!item.equals("-")) q.enqueue(item);
            else if (!q.isEmpty()) StdOut.print(q.dequeue() + " ");
        }
        StdOut.println("(" + q.size() + " left on queue)");
    }   
}
