package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/5/17.
 */
/******************************************************************************
 *  Compilation:  javac EdgeWeightedGraph.java
 *  Execution:    java EdgeWeightedGraph filename.txt
 *  Dependencies: Bag.java Edge.java In.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/43mst/tinyEWG.txt
 *                http://algs4.cs.princeton.edu/43mst/mediumEWG.txt
 *                http://algs4.cs.princeton.edu/43mst/largeEWG.txt
 *
 *  An edge-weighted undirected graph, implemented using adjacency lists.
 *  Parallel edges and self-loops are permitted.
 *
 *  % java EdgeWeightedGraph tinyEWG.txt
 *  8 16
 *  0: 6-0 0.58000  0-2 0.26000  0-4 0.38000  0-7 0.16000
 *  1: 1-3 0.29000  1-2 0.36000  1-7 0.19000  1-5 0.32000
 *  2: 6-2 0.40000  2-7 0.34000  1-2 0.36000  0-2 0.26000  2-3 0.17000
 *  3: 3-6 0.52000  1-3 0.29000  2-3 0.17000
 *  4: 6-4 0.93000  0-4 0.38000  4-7 0.37000  4-5 0.35000
 *  5: 1-5 0.32000  5-7 0.28000  4-5 0.35000
 *  6: 6-4 0.93000  6-0 0.58000  3-6 0.52000  6-2 0.40000
 *  7: 2-7 0.34000  1-7 0.19000  0-7 0.16000  5-7 0.28000  4-7 0.37000
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.Stack;
import edu.princeton.cs.algs4.StdOut;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 *  The {@code EdgeWeightedGraph} class represents an edge-weighted
 *  graph of vertices named 0 through <em>V</em> - 1, where each
 *  undirected edge is of type {@link Edge} and has a real-valued weight.
 *  It supports the following two primary operations: add an edge to the graph,
 *  iterate over all of the edges incident to a vertex. It also provides
 *  methods for returning the number of vertices <em>V</em> and the number
 *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which
 *  is a vertex-indexed array of @link{Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the edges incident to a given vertex, which takes
 *  time proportional to the number of such edges.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/43mst">Section 4.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class EdgeWeightedGraph {
  private static final String NEWLINE = System.getProperty("line.separator");

  private final int V;
  private int E;
  //prefer List over array when used with generics
  private List<Bag<Edge>> adj;

  /**
   * initialize an edge-weighted graph with specific number of vertices
   * @param V
   */
  public EdgeWeightedGraph(int V) {
    if (V < 0) throw new IllegalArgumentException("number of vertices must be nonnegative");
    this.V = V;
    this.E = 0;
    adj = new ArrayList<>(V);
    for (int i = 0; i < V; i++)
      adj.add(new Bag<>());
  }
  /**
   * Initializes a new edge-weighted graph that is a deep copy of {@code G}.
   *
   * @param  G the edge-weighted graph to copy
   */
  public EdgeWeightedGraph(EdgeWeightedGraph G) {
    this(G.V());
    this.E = G.E();
    for (int v = 0; v < G.V(); v++) {
      // reverse so that adjacency list is in same order as original
      Stack<Edge> reverse = new Stack<Edge>();
      for (Edge e : G.adj.get(v)) {
        reverse.push(e);
      }
      for (Edge e : reverse) {
        adj.get(v).add(e);
      }
    }
  }
  /**
   * Initializes an edge-weighted graph from an input stream.
   * The format is the number of vertices <em>V</em>,
   * followed by the number of edges <em>E</em>,
   * followed by <em>E</em> pairs of vertices and edge weights,
   * with each entry separated by whitespace.
   *
   * @param  in the input stream
   * @throws IllegalArgumentException if the endpoints of any edge are not in prescribed range
   * @throws IllegalArgumentException if the number of vertices or edges is negative
   */
  public EdgeWeightedGraph(In in) {
    this(in.readInt());
    int E = in.readInt();
    if (E < 0) throw new IllegalArgumentException("Number of edges must be nonnegative");
    for (int i = 0; i < E; i++) {
      int v = in.readInt();
      int w = in.readInt();
      validateVertex(v);
      validateVertex(w);
      double weight = in.readDouble();
      Edge e = new Edge(v, w, weight);
      addEdge(e);
    }
  }
  public int V() {
    return V;
  }
  public int E() {
    return E;
  }

  /**
   *
   * @param v
   * @return
   */
  public Iterable<Edge> adj(int v) {
    //it is okay to return the original object
    //because iterable has no method to
    //change itself
    return adj.get(v);
  }
  private void addEdge(Edge e) {
    int s = e.either();
    int t = e.other(s);
    adj.get(s).add(e);
    adj.get(t).add(e);
    E++;
  }
  /**
   * Returns the degree of vertex {@code v}.
   *
   * @param  v the vertex
   * @return the degree of vertex {@code v}
   * @throws IllegalArgumentException unless {@code 0 <= v < V}
   */
  public int degree(int v) {
    validateVertex(v);
    return adj.get(v).size();
  }
  /**
   * Returns all edges in this edge-weighted graph.
   * To iterate over the edges in this edge-weighted graph, use foreach notation:
   * {@code for (Edge e : G.edges())}.
   *
   * for performance purposes, we should have a private class
   * implementing Iterable interface, but for now we ignore it
   *
   * @return all edges in this edge-weighted graph, as an iterable
   */
  public Iterable<Edge> edges() {
    List<Edge> allEdges = new ArrayList<>(E);
    int selfLoops = 0;
    for (int i = 0; i < V; i++) {
      for (Edge e: adj.get(i)) {
        if (i < e.other(i)) {
          allEdges.add(e);
        } else if (i == e.other(i)) {
          selfLoops++;
          if (selfLoops % 2 == 0)
            allEdges.add(e);
        }
      }
    }
    return allEdges;
  }

  /**
   * Returns a string representation of the edge-weighted graph.
   * This method takes time proportional to <em>E</em> + <em>V</em>.
   *
   * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
   *         followed by the <em>V</em> adjacency lists of edges
   */
  public String toString() {
    StringBuilder s = new StringBuilder();
    s.append(V + " " + E + NEWLINE);
    for (int v = 0; v < V; v++) {
      s.append(v + ": ");
      for (Edge e : adj.get(v)) {
        s.append(e + "  ");
      }
      s.append(NEWLINE);
    }
    return s.toString();
  }

  // throw an IllegalArgumentException unless {@code 0 <= v < V}
  private void validateVertex(int v) {
    if (v < 0 || v >= V)
      throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
  }
  /**
   * Unit tests the {@code EdgeWeightedGraph} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    EdgeWeightedGraph G = new EdgeWeightedGraph(in);
    StdOut.println(G);
  }
}
