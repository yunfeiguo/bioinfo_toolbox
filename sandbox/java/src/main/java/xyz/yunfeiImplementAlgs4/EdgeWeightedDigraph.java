package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/24/17.
 */
/******************************************************************************
 *  Compilation:  javac EdgeWeightedDigraph.java
 *  Execution:    java EdgeWeightedDigraph digraph.txt
 *  Dependencies: Bag.java DirectedEdge.java
 *  Data files:   http://algs4.cs.princeton.edu/44st/tinyEWD.txt
 *                http://algs4.cs.princeton.edu/44st/mediumEWD.txt
 *                http://algs4.cs.princeton.edu/44st/largeEWD.txt
 *
 *  An edge-weighted digraph, implemented using adjacency lists.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

import java.util.ArrayList;
import java.util.List;

/**
 *  The {@code EdgeWeightedDigraph} class represents a edge-weighted
 *  digraph of vertices named 0 through <em>V</em> - 1, where each
 *  directed edge is of type {@link DirectedEdge} and has a real-valued weight.
 *  It supports the following two primary operations: add a directed edge
 *  to the digraph and iterate over all of edges incident from a given vertex.
 *  It also provides
 *  methods for returning the number of vertices <em>V</em> and the number
 *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
 *  <p>
 *  This implementation uses an adjacency-lists representation, which
 *  is a vertex-indexed array of @link{Bag} objects.
 *  All operations take constant time (in the worst case) except
 *  iterating over the edges incident from a given vertex, which takes
 *  time proportional to the number of such edges.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/44sp">Section 4.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class EdgeWeightedDigraph {
  private static final String NEWLINE = "\n";
  private final int V;
  private int E;
  private List<List<DirectedEdge>> adj = new ArrayList<>();
  private List<Integer> inDegree = new ArrayList<>();

  /**
   * initialize empty graph with V nodes
   * @param V
   * @throws IllegalArgumentException if V < 0
   */
  public EdgeWeightedDigraph(int V) {
    if (V < 0) throw new IllegalArgumentException();
    this.V = V;
    E = 0;
    for (int i = 0; i < V; i++) {
      adj.add(new ArrayList<DirectedEdge>());
      inDegree.add(0);
    }
  }

  /**
   * initilize random graph with V nodes and E random edges
   * @param V
   * @param E
   */
  public EdgeWeightedDigraph(int V, int E) {
    this(V);
    for (int i = 0; i < E; i++) {
      int v = StdRandom.uniform(V);
      int w = StdRandom.uniform(V);
      addEdge(new DirectedEdge(v, w, StdRandom.uniform(100) * 0.01));
    }
  }

  /**
   * make a copy of {@code g}
   * @param g
   */
  public EdgeWeightedDigraph(EdgeWeightedDigraph g) {
    this(g.V());
    for (DirectedEdge e : g.edges()) {
      addEdge(e);
    }
  }
  public EdgeWeightedDigraph(In in) {
    this(in.readInt());
    int E = in.readInt();
    if (E < 0) throw new IllegalArgumentException("Number of edges must be nonnegative");
    for (int i = 0; i < E; i++) {
      int v = in.readInt();
      int w = in.readInt();
      validateVertex(v);
      validateVertex(w);
      double weight = in.readDouble();
      addEdge(new DirectedEdge(v, w, weight));
    }
  }
  public Iterable<DirectedEdge> adj(int v) {
    validateVertex(v);
    return adj.get(v);
  }
  public int V() {
    return V;
  }
  public int E() {
    return E;
  }
  public void addEdge(DirectedEdge e) {
    int from = e.from();
    int to = e.to();
    validateVertex(from);
    validateVertex(to);
    adj.get(from).add(e);
    inDegree.set(to, inDegree.get(to) + 1);
    E++;
  }
  public int outDegree(int v) {
    validateVertex(v);
    return adj.get(v).size();
  }
  public int inDegree(int v) {
    validateVertex(v);
    return inDegree.get(v);
  }
  public Iterable<DirectedEdge> edges() {
    List<DirectedEdge> allEdges = new ArrayList<>();
    for (List<DirectedEdge> l : adj) {
      allEdges.addAll(l);
    }
    return allEdges;
  }
  private void validateVertex(int v) {
    if (v < 0 || v >= V) throw new IllegalArgumentException("not in graph");
  }
  /**
   * Returns a string representation of this edge-weighted digraph.
   *
   * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
   *         followed by the <em>V</em> adjacency lists of edges
   */
  @Override
  public String toString() {
    StringBuilder s = new StringBuilder();
    s.append(V + " " + E + NEWLINE);
    for (int v = 0; v < V; v++) {
      s.append(v + ": ");
      for (DirectedEdge e : adj.get(v)) {
        s.append(e + "  ");
      }
      s.append(NEWLINE);
    }
    return s.toString();
  }
  /**
   * Unit tests the {@code EdgeWeightedDigraph} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    EdgeWeightedDigraph G = new EdgeWeightedDigraph(in);
    StdOut.println(G);
  }
}
