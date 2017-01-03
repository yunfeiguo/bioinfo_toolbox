package yunfeiImplementAlgs4;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

/******************************************************************************
 *  Compilation:  javac AdjMatrixGraph.java
 *  Execution:    java AdjMatrixGraph V E
 *  Dependencies: StdOut.java
 *
 *  A graph, implemented using an adjacency matrix.
 *  Parallel edges are disallowed; self-loops are allowd.
 *  
 ******************************************************************************/

public class AdjMatrixGraph {
    private static final String NEWLINE = System.getProperty("line.separator");
    private int V;
    private int E;
    private boolean[][] adj;

    public AdjMatrixGraph(int V) {
        if (V < 0)
            throw new RuntimeException("number of vertices must be positive");
        this.V = V;
        this.E = 0;
        this.adj = new boolean[V][V];
    }
    //random graph with V vertices and E edges
    public AdjMatrixGraph(int V, int E) {
        this(V);
        if (E < 0)
            throw new RuntimeException("number of edges must be positive");
        if (E > V*V) throw new RuntimeException("Too many edges");

        while (this.E != E) {
            int v = StdRandom.uniform(V);
            int w = StdRandom.uniform(V);
            addEdge(v, w);            
        }
    }
    public int V() {
        return V;
    }
    public int E() {
        return E;
    }
    //add undirected edge v-w
    public void addEdge(int v, int w) {
        if (!adj[v][w]) E++;
        adj[v][w] = true;
        adj[w][v] = true;
    }
    public boolean contains(int v, int w) {
        return adj[v][w];
    }
    public Iterable<Integer> adj(int v) {
        return new AdjIterator(v);
    }
    /*
    private class AdjIterator implements Iterable<Integer>{
        private List<Integer> neighbors;
        public AdjIterator(int v) {
            this.neighbors = new ArrayList<Integer>();
            for (int i = 0; i < V; i++) {
                if (adj[v][i]) {
                    neighbors.add(i);
                }
            }            
        }
        public Iterator<Integer> iterator() {
            return neighbors.iterator();
        }

    }*/
    /*alternative iteratable implementation*/
    private class AdjIterator implements Iterator<Integer>, Iterable<Integer> {
        private int v;
        private int w = 0;

        public AdjIterator(int v) {
            this.v = v;
        }
        public Iterator<Integer> iterator() {
            return this;
        }
        public boolean hasNext() {
            while (w < V) {
                if (adj[v][w]) return true;
                w++;
            }
            return false;
        }
        public Integer next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            return w++;
        }
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
    // string representation of Graph - takes quadratic time
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(V + " " + E + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (int w : adj(v)) {
                s.append(w + " ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    } 


    // test client
    public static void main(String[] args) {
        int V = Integer.parseInt(args[0]);
        int E = Integer.parseInt(args[1]);
        AdjMatrixGraph G = new AdjMatrixGraph(V, E);
        StdOut.println(G);
    }
}
