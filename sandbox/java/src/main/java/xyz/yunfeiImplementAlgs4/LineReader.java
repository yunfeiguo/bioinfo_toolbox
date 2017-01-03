package yunfeiImplementAlgs4;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public class LineReader implements Iterable<String>{
    private BufferedReader br;
    public LineReader(String filename) {
        try {
            br = new BufferedReader(new FileReader(filename));            
        } catch(FileNotFoundException e) {
            System.out.println("File not found");
            e.printStackTrace();
        }
    }
    @Override
    public Iterator<String> iterator() {
        return new LineReaderIterator();
    }
    private class LineReaderIterator implements Iterator<String> {
        String nextElement;
        boolean nextCalled; //record whether next() has been called after hasNext()
        public LineReaderIterator() {
            nextElement = null;
            nextCalled = true; //if false, hasNext will not proceed
        }
        public String next() {
            if (hasNext()) {
                nextCalled = true;
                return nextElement;
            } else {
                throw new NoSuchElementException();
            }
        }
        public boolean hasNext() {
            if(nextCalled) {
                try {
                    nextElement = br.readLine();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                nextCalled = false;
            }
            //if no more elements we close stream
            if (nextElement == null) {
                try {
                    br.close();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
            return nextElement != null;
        }
        public void remove() {
            throw new NotImplementedException();
        }
    }    
    public static void main(String[] args) {
        for (String s : new LineReader("routes.txt")) {
            System.out.println(s);
        }
    }
}
