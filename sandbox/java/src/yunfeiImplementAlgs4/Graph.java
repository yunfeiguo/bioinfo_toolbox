package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac Graph.java        
 *  Execution:    java Graph input.txt
 *  Dependencies: Bag.java In.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/41graph/tinyG.txt
 *
 *  A graph, implemented using an array of sets.
 *  Parallel edges and self-loops allowed.
 *
 *  % java Graph tinyG.txt
 *  13 vertices, 13 edges 
 *  0: 6 2 1 5 
 *  1: 0 
 *  2: 0 
 *  3: 5 4 
 *  4: 5 6 3 
 *  5: 3 4 0 
 *  6: 0 4 
 *  7: 8 
 *  8: 7 
 *  9: 11 10 12 
 *  10: 9 
 *  11: 9 12 
 *  12: 11 9 
 *
 *  % java Graph mediumG.txt
 *  250 vertices, 1273 edges 
 *  0: 225 222 211 209 204 202 191 176 163 160 149 114 97 80 68 59 58 49 44 24 15 
 *  1: 220 203 200 194 189 164 150 130 107 72 
 *  2: 141 110 108 86 79 51 42 18 14 
 *  ...
 *  
 ******************************************************************************/
/**
 *  The <tt>Graph</tt> class represents an undirected graph of vertices
 *  named 0 through <em>V</em> - 1.
 *  It supports the following two primary operations: add an edge to the graph,
 *  iterate over all of the vertices adjacent to a vertex. It also provides
 *  methods for returning the number of vertices <em>V</em> and the number
 *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which 
 *  is a vertex-indexed array of {@link Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the vertices adjacent to a given vertex, which takes
 *  time proportional to the number of such vertices.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41graph">Section 4.1</a>
 *  of <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Graph {
    private static final String NEWLINE = System.getProperty("line.separator");
    
    private final int V; //why final here? because our implementation requires a fixed number of vertices
    private int E;
    private Bag<Integer>[] adj;
    /**
     * Initializes an empty graph with <tt>V</tt> vertices and 0 edges
     * @param V number of vertices
     * @throws IllegalArgumentException if <tt>V</tt> < 0
     */
    @SuppressWarnings("unchecked")
    public Graph(int V) {
        if (V < 0)
            throw new IllegalArgumentException();
        this.V = V;
        this.E = 0;
        //adj = new Bag<Integer>[V]; //wrong code, cannot declare a generic type of array
        adj = (Bag<Integer>[]) new Bag[V];
        for (int v = 0; v < V; v++) {
            adj[v] = new Bag<Integer>();
        }        
    }
    /**
     * initializes a graph from an input stream.
     * The format is the number of vertices <em>V</em>,
     * followed by the number of edges <em>E</em>,
     * followed by <em>E</em> pairs of vertices, with each entry separated by whitespace.
     * 
     * @param in input stream
     * @throws IndexOutOfBoundsException
     * @throws {@link IllegalArgumentException} 
     */
    public Graph(In in) {
        this(in.readInt());
        int E = in.readInt();
        if (E < 0) 
            throw new IllegalArgumentException();
        for (int i = 0; i < E; i++) {
            int v = in.readInt();
            int w = in.readInt();
            addEdge(v, w);
        }
    }
    /**
     * initializes a new graph that is a deep copy of <tt>G</tt>
     * @param G the graph to copy
     */
    public Graph(Graph G) {
        this(G.V());
        for (int i = 0; i < G.V(); i++) {
            for (int j : G.adj(i))
                addEdge(i, j);
        }
    }
    /**
     * 
     * @return number of vertices
     */
    public int V() {
        return this.V;
    }
    /**
     * 
     * @return number of edges
     */
    public int E() {
        return this.E;
    }
    public Iterable<Integer> adj(int v) {
        return adj[v];
    }
    public void addEdge(int v, int w) {
        validateVertex(v);
        validateVertex(w);
        E++;
        adj[v].add(w);
        adj[w].add(v);
    }
    private void validateVertex(int v) {
        if (v < 0 || v >= V) 
            throw new IndexOutOfBoundsException();        
    }
    /**
     * 
     * @param v the vertex being queried
     * @return the degree of vertex
     */
    public int degree(int v) {
        validateVertex(v);
        return adj[v].size();
    }
    
    /**
     * @return number of vertices, followed by number of edges E followed by V adjacency lists 
     */
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(V + " vertices, " + E + " edges " + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (int w : adj[v]) {
                s.append(w + " ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }
    /**
     * Unit tests the <tt>Graph</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Graph G = new Graph(in);
        StdOut.println(G);
    }

}
