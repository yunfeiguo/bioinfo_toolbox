package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/25/17.
 */
/******************************************************************************
 *  Compilation:  javac DijkstraSP.java
 *  Execution:    java DijkstraSP input.txt s
 *  Dependencies: EdgeWeightedDigraph.java IndexMinPQ.java Stack.java DirectedEdge.java
 *  Data files:   http://algs4.cs.princeton.edu/44sp/tinyEWD.txt
 *                http://algs4.cs.princeton.edu/44sp/mediumEWD.txt
 *                http://algs4.cs.princeton.edu/44sp/largeEWD.txt
 *
 *  Dijkstra's algorithm. Computes the shortest path tree.
 *  Assumes all weights are nonnegative.
 *
 *  % java DijkstraSP tinyEWD.txt 0
 *  0 to 0 (0.00)
 *  0 to 1 (1.05)  0->4  0.38   4->5  0.35   5->1  0.32
 *  0 to 2 (0.26)  0->2  0.26
 *  0 to 3 (0.99)  0->2  0.26   2->7  0.34   7->3  0.39
 *  0 to 4 (0.38)  0->4  0.38
 *  0 to 5 (0.73)  0->4  0.38   4->5  0.35
 *  0 to 6 (1.51)  0->2  0.26   2->7  0.34   7->3  0.39   3->6  0.52
 *  0 to 7 (0.60)  0->2  0.26   2->7  0.34
 *
 *  % java DijkstraSP mediumEWD.txt 0
 *  0 to 0 (0.00)
 *  0 to 1 (0.71)  0->44  0.06   44->93  0.07   ...  107->1  0.07
 *  0 to 2 (0.65)  0->44  0.06   44->231  0.10  ...  42->2  0.11
 *  0 to 3 (0.46)  0->97  0.08   97->248  0.09  ...  45->3  0.12
 *  0 to 4 (0.42)  0->44  0.06   44->93  0.07   ...  77->4  0.11
 *  ...
 *
 ******************************************************************************/


import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import org.omg.PortableInterceptor.INACTIVE;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

/**
 *  The {@code DijkstraSP} class represents a data type for solving the
 *  single-source shortest paths problem in edge-weighted digraphs
 *  where the edge weights are nonnegative.
 *  <p>
 *  This implementation uses Dijkstra's algorithm with a binary heap.
 *  The constructor takes time proportional to <em>E</em> log <em>V</em>,
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the {@code distTo()} and {@code hasPathTo()} methods take
 *  constant time and the {@code pathTo()} method takes time proportional to the
 *  number of edges in the shortest path returned.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/44sp">Section 4.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class DijkstraSP {
  private List<Double> distanceTo;
  private List<DirectedEdge> directParentOf;
  //used in intermediate calculation
  private IndexMinPQ<Double> dynamicDistanceTo;
  private int source;
  /**
   * initialize instance and calculate shortest paths
   * from source to all other nodes
   * @param g
   * @param source
   */
  public DijkstraSP(EdgeWeightedDigraph g, int source) {
    this.source = source;
    distanceTo = new ArrayList<>(g.V());
    directParentOf = new ArrayList<>(g.V());
    dynamicDistanceTo = new IndexMinPQ<>(g.V());
    for (int i = 0; i < g.V(); i++) {
      distanceTo.add(Double.NaN);
      directParentOf.add(null);
      dynamicDistanceTo.insert(i, Double.POSITIVE_INFINITY);
    }
    validateVertex(source);
    //set source distance to 0
    dynamicDistanceTo.decreaseKey(source, (double) 0);
    directParentOf.set(source, null);
    distanceTo.set(source, (double) 0);

    while(!dynamicDistanceTo.isEmpty()) {
      int next = dynamicDistanceTo.delMin();
      for (DirectedEdge e : g.adj(next)) {
        relax(e);
      }
    }

    // check optimality conditions
    assert check(g, source);
  }

  /**
   * return shortest path from source to destination
   * @param destination
   * @return
   */
  public Iterable<DirectedEdge> pathTo(int destination) {
    if (!hasPathTo(destination)) return new ArrayList<>();
    Deque<DirectedEdge> path = new ArrayDeque<>();
    /*
    in implementation, it is better to determine whether edge is null than to determine whether
    we have reached source vertex
     */
    for(DirectedEdge e = directParentOf.get(destination); e != null; e = directParentOf.get(e.from())) {
        path.addFirst(e);
    }
    return path;
  }
  public boolean hasPathTo(int destination) {
    validateVertex(destination);
    return !Double.isNaN(distanceTo.get(destination));
  }
  public double distTo(int destination) {
    if (!hasPathTo(destination)) throw new IllegalArgumentException("no path");
    return distanceTo.get(destination);
  }
  private void relax(DirectedEdge e) {
    int from = e.from();
    int to = e.to();
    if (Double.isNaN(distanceTo.get(to)) || distanceTo.get(from) + e.weight() < distanceTo.get(to)) {
      distanceTo.set(to, distanceTo.get(from) + e.weight());
      directParentOf.set(to, e);
      dynamicDistanceTo.decreaseKey(to, distanceTo.get(to));
    }
  }
  // check optimality conditions:
  // (i) for all edges e:            distTo[e.to()] <= distTo[e.from()] + e.weight()
  // (ii) for all edge e on the SPT: distTo[e.to()] == distTo[e.from()] + e.weight()
  private boolean check(EdgeWeightedDigraph G, int s) {

    // check that edge weights are nonnegative
    for (DirectedEdge e : G.edges()) {
      if (e.weight() < 0) {
        System.err.println("negative edge weight detected");
        return false;
      }
    }

    // check that distTo[v] and edgeTo[v] are consistent
    if (!Double.isNaN(distanceTo.get(s)) || directParentOf.get(s) != null) {
      System.err.println("distTo[s] and edgeTo[s] inconsistent");
      return false;
    }
    for (int v = 0; v < G.V(); v++) {
      if (v == s) continue;
      if (directParentOf.get(v) == null && distanceTo.get(v) != Double.POSITIVE_INFINITY) {
        System.err.println("distTo[] and edgeTo[] inconsistent");
        return false;
      }
    }

    // check that all edges e = v->w satisfy distTo[w] <= distTo[v] + e.weight()
    for (int v = 0; v < G.V(); v++) {
      for (DirectedEdge e : G.adj(v)) {
        int w = e.to();
        if (distanceTo.get(v) + e.weight() < distanceTo.get(w)) {
          System.err.println("edge " + e + " not relaxed");
          return false;
        }
      }
    }

    // check that all edges e = v->w on SPT satisfy distTo[w] == distTo[v] + e.weight()
    for (int w = 0; w < G.V(); w++) {
      if (directParentOf.get(w) == null) continue;
      DirectedEdge e = directParentOf.get(w);
      int v = e.from();
      if (w != e.to()) return false;
      if (distanceTo.get(w) + e.weight() != distanceTo.get(w)) {
        System.err.println("edge " + e + " on shortest path not tight");
        return false;
      }
    }
    return true;
  }

  // throw an IllegalArgumentException unless {@code 0 <= v < V}
  private void validateVertex(int v) {
    int V = distanceTo.size();
    if (v < 0 || v >= V)
      throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
  }

  /**
   * Unit tests the {@code DijkstraSP} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    EdgeWeightedDigraph G = new EdgeWeightedDigraph(in);
    int s = Integer.parseInt(args[1]);

    // compute shortest paths
    DijkstraSP sp = new DijkstraSP(G, s);


    // print shortest path
    for (int t = 0; t < G.V(); t++) {
      if (sp.hasPathTo(t)) {
        StdOut.printf("%d to %d (%.2f)  ", s, t, sp.distTo(t));
        for (DirectedEdge e : sp.pathTo(t)) {
          StdOut.print(e + "   ");
        }
        StdOut.println();
      }
      else {
        StdOut.printf("%d to %d         no path\n", s, t);
      }
    }
  }
}
