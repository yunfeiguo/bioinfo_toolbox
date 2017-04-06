package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/27/17.
 */
/******************************************************************************
 *  Compilation:  javac EdgeWeightedDirectedCycle.java
 *  Execution:    java EdgeWeightedDirectedCycle V E F
 *  Dependencies: EdgeWeightedDigraph.java DirectedEdge.java Stack.java
 *
 *  Finds a directed cycle in an edge-weighted digraph.
 *  Runs in O(E + V) time.
 *
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

/**
 *  The {@code EdgeWeightedDirectedCycle} class represents a data type for
 *  determining whether an edge-weighted digraph has a directed cycle.
 *  The <em>hasCycle</em> operation determines whether the edge-weighted
 *  digraph has a directed cycle and, if so, the <em>cycle</em> operation
 *  returns one.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>hasCycle</em> operation takes constant time;
 *  the <em>cycle</em> operation takes time proportional
 *  to the length of the cycle.
 *  <p>
 *  See {@link Topological} to compute a topological order if the edge-weighted
 *  digraph is acyclic.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/44sp">Section 4.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class EdgeWeightedDirectedCycle {
  private List<Boolean> isBeingVisited;
  private List<Boolean> visited;
  private List<DirectedEdge> edgeTo;
  private Deque<DirectedEdge> cycle;
  public EdgeWeightedDirectedCycle(EdgeWeightedDigraph g) {
    isBeingVisited = new ArrayList<>(g.V());
    edgeTo = new ArrayList<>(g.V());
    visited = new ArrayList<>(g.V());
    for (int i = 0; i < g.V(); i++) {
      isBeingVisited.add(false);
      visited.add(false);
      edgeTo.add(null);
    }

    for (int i = 0; i < g.V(); i++) {
      dfs(g, i);
    }
  }
  public boolean hasCycle() {
    return cycle != null;
  }
  public Iterable<DirectedEdge> cycle() {
    return cycle;
  }
  private void dfs(EdgeWeightedDigraph g, int v) {
    isBeingVisited.set(v, true);
    for (DirectedEdge e : g.adj(v)) {
      //why put condition here? because we could find a cycle
      //before iterating over all edges of a parent vertex
      //we must avoid any more traversal after finding a cycle
      if (hasCycle() || visited.get(e.to())) return;
      int to = e.to();
      if (isBeingVisited.get(to)) {
        cycle = new ArrayDeque<>();
        for (; e.from() != to; e = edgeTo.get(e.from())) {
          cycle.addFirst(e);
        }
        cycle.addFirst(e);
        return;
      }
      if (!visited.get(to)) {
        edgeTo.set(to, e);
        dfs(g, to);
      }
    }
    isBeingVisited.set(v, false);
    visited.set(v, true);
  }
  /**
   * Unit tests the {@code EdgeWeightedDirectedCycle} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {

    // create random DAG with V vertices and E edges; then add F random edges
    int V = 3;
    int E = 3;
    int F = 3;
    StdRandom.setSeed(11l);
    EdgeWeightedDigraph G = new EdgeWeightedDigraph(V);
    int[] vertices = new int[V];
    for (int i = 0; i < V; i++)
      vertices[i] = i;
    StdRandom.shuffle(vertices);
    for (int i = 0; i < E; i++) {
      int v, w;
      do {
        v = StdRandom.uniform(V);
        w = StdRandom.uniform(V);
      } while (v >= w);
      double weight = StdRandom.uniform();
      G.addEdge(new DirectedEdge(v, w, weight));
    }

    // add F extra edges
    for (int i = 0; i < F; i++) {
      int v = StdRandom.uniform(V);
      int w = StdRandom.uniform(V);
      double weight = StdRandom.uniform(0.0, 1.0);
      G.addEdge(new DirectedEdge(v, w, weight));
    }

    StdOut.println(G);

    // find a directed cycle
    EdgeWeightedDirectedCycle finder = new EdgeWeightedDirectedCycle(G);
    if (finder.hasCycle()) {
      StdOut.print("Cycle: ");
      for (DirectedEdge e : finder.cycle()) {
        StdOut.print(e + " ");
      }
      StdOut.println();
    }

    // or give topologial sort
    else {
      StdOut.println("No directed cycle");
    }
  }
}
