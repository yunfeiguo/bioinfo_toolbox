package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/27/17.
 */
/******************************************************************************
 *  Compilation:  javac AcyclicSP.java
 *  Execution:    java AcyclicSP V E
 *  Dependencies: EdgeWeightedDigraph.java DirectedEdge.java Topological.java
 *  Data files:   http://algs4.cs.princeton.edu/44sp/tinyEWDAG.txt
 *
 *  Computes shortest paths in an edge-weighted acyclic digraph.
 *
 *  % java AcyclicSP tinyEWDAG.txt 5
 *  5 to 0 (0.73)  5->4  0.35   4->0  0.38
 *  5 to 1 (0.32)  5->1  0.32
 *  5 to 2 (0.62)  5->7  0.28   7->2  0.34
 *  5 to 3 (0.61)  5->1  0.32   1->3  0.29
 *  5 to 4 (0.35)  5->4  0.35
 *  5 to 5 (0.00)
 *  5 to 6 (1.13)  5->1  0.32   1->3  0.29   3->6  0.52
 *  5 to 7 (0.28)  5->7  0.28
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

/**
 *  The {@code AcyclicSP} class represents a data type for solving the
 *  single-source shortest paths problem in edge-weighted directed acyclic
 *  graphs (DAGs). The edge weights can be positive, negative, or zero.
 *  <p>
 *  This implementation uses a topological-sort based algorithm.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>,
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
public class AcyclicSP {
  private final List<Double> distTo;
  private final List<DirectedEdge> edgeTo;

  public AcyclicSP(EdgeWeightedDigraph g, int source) {
    distTo = new ArrayList<>(g.V());
    edgeTo = new ArrayList<>(g.V());

    for (int i = 0; i < g.V(); i++) {
      distTo.add(Double.POSITIVE_INFINITY);
      edgeTo.add(null);
    }
    //if we do not set source distance to 0
    //relax() will have no effect
    distTo.set(source, (double) 0);

    Topological topoSortedDigraph = new Topological(g);
    for (int i : topoSortedDigraph.order()) {
      for (DirectedEdge e : g.adj(i)) {
        relax(e);
      }
    }
  }

  /**
   * update distTo and edgeTo if we find a shorter path
   * via e
   * @param e
   */
  private void relax(DirectedEdge e) {
    int to = e.to();
    int from = e.from();
    if (distTo.get(to) > distTo(from) + e.weight()) {
      edgeTo.set(to, e);
      distTo.set(to, distTo.get(from) + e.weight());
    }
  }
  public double distTo(int v) {
    validateVertex(v);
    return distTo.get(v);
  }
  // throw an IllegalArgumentException unless {@code 0 <= v < V}
  private void validateVertex(int v) {
    int V = distTo.size();
    if (v < 0 || v >= V)
      throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
  }
  public boolean hasPathTo(int v) {
    validateVertex(v);
    //cannot use this for determination, why? because
    //for source itself, there is no edge to it
    //return edgeTo.get(v) != null;
    return distTo.get(v) < Double.POSITIVE_INFINITY;
  }
  public Iterable<DirectedEdge> pathTo(int v) {
    if (!hasPathTo(v)) throw new IllegalArgumentException("no path");
    Deque<DirectedEdge> path = new ArrayDeque<>();
    for (DirectedEdge e = edgeTo.get(v); e != null; e = edgeTo.get(e.from())) {
      path.addFirst(e);
    }
    return path;
  }

  /**
   * Unit tests the {@code AcyclicSP} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    int s = Integer.parseInt(args[1]);
    EdgeWeightedDigraph G = new EdgeWeightedDigraph(in);

    // find shortest path from s to each other vertex in DAG
    AcyclicSP sp = new AcyclicSP(G, s);
    for (int v = 0; v < G.V(); v++) {
      if (sp.hasPathTo(v)) {
        StdOut.printf("%d to %d (%.2f)  ", s, v, sp.distTo(v));
        for (DirectedEdge e : sp.pathTo(v)) {
          StdOut.print(e + "   ");
        }
        StdOut.println();
      }
      else {
        StdOut.printf("%d to %d         no path\n", s, v);
      }
    }
  }
}
