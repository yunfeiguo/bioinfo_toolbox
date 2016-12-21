package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import java.nio.channels.UnsupportedAddressTypeException;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class FixedCapacityStackOfStrings implements Iterable<String> {
    private String[] a;
    private int N;
    
    public FixedCapacityStackOfStrings(int capacity) {
        a = new String[capacity];
        N = 0;
    }
    
    public void push(String x) {
        if (isFull()) throw new RuntimeException();
        a[N++] = x;
    }
    
    public String pop() {
        if (isEmpty()) throw new NoSuchElementException();
        return a[--N];
    }
    
    public String peek() {
        if (isEmpty()) throw new NoSuchElementException();
        return a[N - 1];
    }
    
    public boolean isEmpty() {
        return N == 0;
    }
    
    public boolean isFull() {
        return N == a.length;
    }
    
    public Iterator<String> iterator() {
        return new ReverseArrayIterator();
    }
    
    public class ReverseArrayIterator implements Iterator<String> {
        private int i = N - 1;
        
        public boolean hasNext() {
            return i >= 0;
        }
        
        public String next() {
            if (!hasNext()) throw new NoSuchElementException();
            return a[i--];
        }
        
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    public static void main(String[] args) {
        int max = Integer.parseInt(args[0]);
        FixedCapacityStackOfStrings stack = new FixedCapacityStackOfStrings(max);
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
