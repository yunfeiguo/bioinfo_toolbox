/******************************************************************************
 *  Compilation:  javac ResizingArrayQueue.java
 *  Execution:    java ResizingArrayQueue < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/13stacks/tobe.txt  
 *  
 *  Queue implementation with a resizing array.
 *
 *  % java ResizingArrayQueue < tobe.txt 
 *  to be or not to be (2 left on queue)
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Stack;

public class ResizingArrayQueue<Item> implements Iterable<Item> {


    /**
     *  The <tt>ResizingArrayQueue</tt> class represents a first-in-first-out (FIFO)
     *  queue of generic items.
     *  It supports the usual <em>enqueue</em> and <em>dequeue</em>
     *  operations, along with methods for peeking at the first item,
     *  testing if the queue is empty, and iterating through
     *  the items in FIFO order.
     *  <p>
     *  This implementation uses a resizing array, which double the underlying array
     *  when it is full and halves the underlying array when it is one-quarter full.
     *  The <em>enqueue</em> and <em>dequeue</em> operations take constant amortized time.
     *  The <em>size</em>, <em>peek</em>, and <em>is-empty</em> operations takes
     *  constant time in the worst case. 
     *  <p>
     *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/13stacks">Section 1.3</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *
     *  @author Robert Sedgewick
     *  @author Kevin Wayne
     */
    private Item[] a;
    private int N;
    private int head;
    private int tail;

    @SuppressWarnings("unchecked")
    public ResizingArrayQueue() {               
        a = (Item[]) new Object[2]; //default capacity is 2
        N = 0;
        head = 0;
        tail = 0;
    }

    public void enqueue(Item x) {
        if (isFull()) {
            resize(a.length * 2);
        }
        a[tail % a.length] = x;
        ++tail;
        ++N;
    }
    @SuppressWarnings("unchecked")
    private void resize(int size) {
        Item[] b = (Item[]) new Object[size];
        //do not use head and tail to bound the index for b
        //when sizing down, index may go out of bound
        for (int i = 0; i < N; i++)
            b[i] = a[(i + head) % a.length];            
        a = b;
        head = 0;
        tail = N;
    }

    public Item dequeue() {
        if (isEmpty()) throw new NoSuchElementException("empty queue");
        Item ret = a[head % a.length];
        a[head % a.length] = null; //avoid loitering, very important!!!
        head++;
        N--;
        if (N > 0 && N <= a.length / 4) resize(a.length / 2);
        return ret;
    }

    public Item peek() {
        if (isEmpty()) throw new NoSuchElementException();
        return a[head % a.length];
    }

    public boolean isEmpty() {
        return N == 0;
    }

    public boolean isFull() {
        return N == a.length;
    }

    public Iterator<Item> iterator() {
        return new ReverseArrayIterator();
    }
    
    public int size() {
        return N;
    }

    public class ReverseArrayIterator implements Iterator<Item> {
        private int i = head;                

        public boolean hasNext() {
            return i < tail;
        }

        public Item next() {
            if (!hasNext()) throw new NoSuchElementException();
            return a[i++ % a.length];
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    public static void main(String[] args) {
        ResizingArrayQueue<String> q = new ResizingArrayQueue<String>();
        while (!StdIn.isEmpty()) {
            String item = StdIn.readString();
            if (!item.equals("-")) q.enqueue(item); 
            else if (q.isEmpty())  StdOut.println("BAD INPUT"); 
            else                       StdOut.print(q.dequeue() + " ");
        }
        StdOut.println();

        // print what's left on the stack
        StdOut.print("Left on queue: ");
        for (String s : q) {
            StdOut.print(s + " ");
        }
        StdOut.println();
    } 
}

