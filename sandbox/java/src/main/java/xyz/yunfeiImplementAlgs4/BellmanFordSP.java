package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/2/17.
 */
/******************************************************************************
 *  Compilation:  javac BellmanFordSP.java
 *  Execution:    java BellmanFordSP filename.txt s
 *  Dependencies: EdgeWeightedDigraph.java DirectedEdge.java Queue.java
 *                EdgeWeightedDirectedCycle.java
 *  Data files:   http://algs4.cs.princeton.edu/44sp/tinyEWDn.txt
 *                http://algs4.cs.princeton.edu/44sp/mediumEWDnc.txt
 *
 *  Bellman-Ford shortest path algorithm. Computes the shortest path tree in
 *  edge-weighted digraph G from vertex s, or finds a negative cost cycle
 *  reachable from s.
 *
 *  % java BellmanFordSP tinyEWDn.txt 0
 *  0 to 0 ( 0.00)
 *  0 to 1 ( 0.93)  0->2  0.26   2->7  0.34   7->3  0.39   3->6  0.52   6->4 -1.25   4->5  0.35   5->1  0.32
 *  0 to 2 ( 0.26)  0->2  0.26
 *  0 to 3 ( 0.99)  0->2  0.26   2->7  0.34   7->3  0.39
 *  0 to 4 ( 0.26)  0->2  0.26   2->7  0.34   7->3  0.39   3->6  0.52   6->4 -1.25
 *  0 to 5 ( 0.61)  0->2  0.26   2->7  0.34   7->3  0.39   3->6  0.52   6->4 -1.25   4->5  0.35
 *  0 to 6 ( 1.51)  0->2  0.26   2->7  0.34   7->3  0.39   3->6  0.52
 *  0 to 7 ( 0.60)  0->2  0.26   2->7  0.34
 *
 *  % java BellmanFordSP tinyEWDnc.txt 0
 *  4->5  0.35
 *  5->4 -0.66
 *
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

import java.util.*;

/**
 *  The {@code BellmanFordSP} class represents a data type for solving the
 *  single-source shortest paths problem in edge-weighted digraphs with
 *  no negative cycles.
 *  The edge weights can be positive, negative, or zero.
 *  This class finds either a shortest path from the source vertex <em>s</em>
 *  to every other vertex or a negative cycle reachable from the source vertex.
 *  <p>
 *  This implementation uses the Bellman-Ford-Moore algorithm.
 *  The constructor takes time proportional to <em>V</em> (<em>V</em> + <em>E</em>)
 *  in the worst case, where <em>V</em> is the number of vertices and <em>E</em>
 *  is the number of edges.
 *  Afterwards, the {@code distTo()}, {@code hasPathTo()}, and {@code hasNegativeCycle()}
 *  methods take constant time; the {@code pathTo()} and {@code negativeCycle()}
 *  method takes time proportional to the number of edges returned.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/44sp">Section 4.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class BellmanFordSP {
  private List<Double> distTo;
  private List<DirectedEdge> edgeTo;
  private Iterable<DirectedEdge> negativeCycle;
  private int passCount;

  public BellmanFordSP(EdgeWeightedDigraph g, int source) {
    distTo = new ArrayList<>(g.V());
    edgeTo = new ArrayList<>(g.V());
    for (int i = 0; i < g.V(); i++) {
      distTo.add(Double.POSITIVE_INFINITY);
      edgeTo.add(null);
    }
    java.util.Queue<Integer> q = new LinkedList<>();
    distTo.set(source, (double) 0);
    q.offer(source);
    passCount = 0;
    while(!q.isEmpty()) {
      int n = q.size();
      for (int i = 0; i < n; i++) {
        int next = q.poll();
        for (DirectedEdge e : g.adj(next)) {
          if (relax(e)) {
            q.offer(e.to());
          }
        }
      }
      passCount++;
      if (passCount % g.V() == 0) {
        //there must exist a negative weight cycle
        //in the edges we chose for shortest path
        EdgeWeightedDigraph subGWithCycle = new EdgeWeightedDigraph(g.V());
        for (DirectedEdge e : edgeTo) {
          if (e != null) {
            subGWithCycle.addEdge(e);
          }
        }
        EdgeWeightedDirectedCycle cycleFinder = new EdgeWeightedDirectedCycle(subGWithCycle);
        negativeCycle = cycleFinder.cycle();
        break;
      }
    }
  }

  /**
   *
   * @param e
   * @return true if update occurs
   */
  private boolean relax(DirectedEdge e) {
    int from = e.from();
    int to = e.to();
    if (distTo.get(to) > distTo.get(from) + e.weight()) {
      distTo.set(to, distTo.get(from) + e.weight());
      edgeTo.set(to, e);
      return true;
    }
    return false;
  }
  private boolean hasNegativeCycle() {
    return negativeCycle != null;
  }
  private Iterable<DirectedEdge> negativeCycle() {
    return negativeCycle;
  }
  public double distTo(int destination) {
    if(!hasPathTo(destination)) throw new IllegalArgumentException("no path");
    return distTo.get(destination);
  }
  public boolean hasPathTo(int destination) {
    return distTo.get(destination) < Double.POSITIVE_INFINITY;
  }
  public Iterable<DirectedEdge> pathTo(int destination) {
    if(!hasPathTo(destination)) throw new IllegalArgumentException("no path");
    Deque<DirectedEdge> path = new ArrayDeque<>();
    for (DirectedEdge e = edgeTo.get(destination); e != null; e = edgeTo.get(e.from())) {
      path.addFirst(e);
    }
    return path;
  }
  // check optimality conditions: either
  // (i) there exists a negative cycle reacheable from s
  //     or
  // (ii)  for all edges e = v->w:            distTo[w] <= distTo[v] + e.weight()
  // (ii') for all edges e = v->w on the SPT: distTo[w] == distTo[v] + e.weight()
  private boolean check(EdgeWeightedDigraph G, int s) {

    // has a negative cycle
    if (hasNegativeCycle()) {
      double weight = 0.0;
      for (DirectedEdge e : negativeCycle()) {
        weight += e.weight();
      }
      if (weight >= 0.0) {
        System.err.println("error: weight of negative cycle = " + weight);
        return false;
      }
    }

    // no negative cycle reachable from source
    else {

      // check that distTo[v] and edgeTo[v] are consistent
      if (distTo.get(s) != 0.0 || edgeTo.get(s) != null) {
        System.err.println("distanceTo[s] and edgeTo[s] inconsistent");
        return false;
      }
      for (int v = 0; v < G.V(); v++) {
        if (v == s) continue;
        if (edgeTo.get(v) == null && distTo.get(v) != Double.POSITIVE_INFINITY) {
          System.err.println("distTo[] and edgeTo[] inconsistent");
          return false;
        }
      }

      // check that all edges e = v->w satisfy distTo[w] <= distTo[v] + e.weight()
      for (int v = 0; v < G.V(); v++) {
        for (DirectedEdge e : G.adj(v)) {
          int w = e.to();
          if (distTo.get(v) + e.weight() < distTo.get(w)) {
            System.err.println("edge " + e + " not relaxed");
            return false;
          }
        }
      }

      // check that all edges e = v->w on SPT satisfy distTo[w] == distTo[v] + e.weight()
      for (int w = 0; w < G.V(); w++) {
        if (edgeTo.get(w) == null) continue;
        DirectedEdge e = edgeTo.get(w);
        int v = e.from();
        if (w != e.to()) return false;
        if (distTo.get(v) + e.weight() != distTo.get(w)) {
          System.err.println("edge " + e + " on shortest path not tight");
          return false;
        }
      }
    }

    StdOut.println("Satisfies optimality conditions");
    StdOut.println();
    return true;
  }

  // throw an IllegalArgumentException unless {@code 0 <= v < V}
  private void validateVertex(int v) {
    int V = distTo.size();
    if (v < 0 || v >= V)
      throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
  }

  /**
   * Unit tests the {@code BellmanFordSP} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    int s = Integer.parseInt(args[1]);
    EdgeWeightedDigraph G = new EdgeWeightedDigraph(in);

    BellmanFordSP sp = new BellmanFordSP(G, s);

    // print negative cycle
    if (sp.hasNegativeCycle()) {
      StdOut.println("has negative cycle");
      for (DirectedEdge e : sp.negativeCycle())
        StdOut.println(e);
    }

    // print shortest paths
    else {
      for (int v = 0; v < G.V(); v++) {
        if (sp.hasPathTo(v)) {
          StdOut.printf("%d to %d (%5.2f)  ", s, v, sp.distTo(v));
          for (DirectedEdge e : sp.pathTo(v)) {
            StdOut.print(e + "   ");
          }
          StdOut.println();
        }
        else {
          StdOut.printf("%d to %d           no path\n", s, v);
        }
      }
    }

  }
}
