/******************************************************************************
 *  Compilation:  javac ResizingArrayStack.java
 *  Execution:    java ResizingArrayStack < input.txt
 *  Dependencies: StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/13stacks/tobe.txt
 *  
 *  Stack implementation with a resizing array.
 *
 *  % more tobe.txt 
 *  to be or not to - be - - that - - - is
 *
 *  % java ResizingArrayStack < tobe.txt
 *  to be not that or be (2 left on stack)
 *
 ******************************************************************************/

package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Stack;

/**
 *  The <tt>ResizingArrayStack</tt> class represents a last-in-first-out (LIFO) stack
 *  of generic items.
 *  It supports the usual <em>push</em> and <em>pop</em> operations, along with methods
 *  for peeking at the top item, testing if the stack is empty, and iterating through
 *  the items in LIFO order.
 *  <p>
 *  This implementation uses a resizing array, which double the underlying array
 *  when it is full and halves the underlying array when it is one-quarter full.
 *  The <em>push</em> and <em>pop</em> operations take constant amortized time.
 *  The <em>size</em>, <em>peek</em>, and <em>is-empty</em> operations takes
 *  constant time in the worst case. 
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/13stacks">Section 1.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class ResizingArrayStack<Item> implements Iterable<Item> {   
        private Item[] a;
        private int N;
        @SuppressWarnings("unchecked")
        public ResizingArrayStack(int capacity) {
            a = (Item[]) new Object[capacity]; //no generic array creation
            N = 0;
            Stack<Integer> s = new Stack<Integer>();
        }
        
        @SuppressWarnings("unchecked")
        public ResizingArrayStack() {
            
            a = (Item[]) new Object[2]; //default capacity is 2
            N = 0;
        }

        public void push(Item x) {
            if (isFull()) {
                resize(a.length * 2);
            }
            a[N++] = x;
        }
        @SuppressWarnings("unchecked")
        private void resize(int size) {
            Item[] b = (Item[]) new Object[size];
            for (int i = 0; i < N; i++)
                b[i] = a[i];            
            a = b;
        }

        public Item pop() {
            if (isEmpty()) throw new NoSuchElementException("stack underflow");
            Item ret = a[--N];
            a[N + 1] = null; //avoid loitering, very important!!!
            if (N <= a.length/4) resize(a.length/2);
            return ret;
        }

        public Item peek() {
            if (isEmpty()) throw new NoSuchElementException();
            return a[N - 1];
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

        public class ReverseArrayIterator implements Iterator<Item> {
            private int i = N - 1;

            public boolean hasNext() {
                return i >= 0;
            }

            public Item next() {
                if (!hasNext()) throw new NoSuchElementException();
                return a[i--];
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        }

        public static void main(String[] args) {
            ResizingArrayStack<String> stack = new ResizingArrayStack<String>();
            while (!StdIn.isEmpty()) {
                String item = StdIn.readString();
                if (!item.equals("-")) stack.push(item); 
                else if (stack.isEmpty())  StdOut.println("BAD INPUT"); 
                else                       StdOut.print(stack.pop() + " ");
            }
            StdOut.println();

            // print what's left on the stack
            StdOut.print("Left on stack: ");
            for (String s : stack) {
                StdOut.print(s + " ");
            }
            StdOut.println();
        } 
    }


