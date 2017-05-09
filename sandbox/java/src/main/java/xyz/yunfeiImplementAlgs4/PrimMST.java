package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/6/17.
 */
/******************************************************************************
 *  Compilation:  javac PrimMST.java
 *  Execution:    java PrimMST filename.txt
 *  Dependencies: EdgeWeightedGraph.java Edge.java Queue.java
 *                IndexMinPQ.java UF.java In.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/43mst/tinyEWG.txt
 *                http://algs4.cs.princeton.edu/43mst/mediumEWG.txt
 *                http://algs4.cs.princeton.edu/43mst/largeEWG.txt
 *
 *  Compute a minimum spanning forest using Prim's algorithm.
 *
 *  %  java PrimMST tinyEWG.txt
 *  1-7 0.19000
 *  0-2 0.26000
 *  2-3 0.17000
 *  4-5 0.35000
 *  5-7 0.28000
 *  6-2 0.40000
 *  0-7 0.16000
 *  1.81000
 *
 *  % java PrimMST mediumEWG.txt
 *  1-72   0.06506
 *  2-86   0.05980
 *  3-67   0.09725
 *  4-55   0.06425
 *  5-102  0.03834
 *  6-129  0.05363
 *  7-157  0.00516
 *  ...
 *  10.46351
 *
 *  % java PrimMST largeEWG.txt
 *  ...
 *  647.66307
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.UF;
import sun.security.x509.EDIPartyName;

import java.util.ArrayList;
import java.util.List;

/**
 *  The {@code PrimMST} class represents a data type for computing a
 *  <em>minimum spanning tree</em> in an edge-weighted graph.
 *  The edge weights can be positive, zero, or negative and need not
 *  be distinct. If the graph is not connected, it computes a <em>minimum
 *  spanning forest</em>, which is the union of minimum spanning trees
 *  in each connected component. The {@code weight()} method returns the
 *  weight of a minimum spanning tree and the {@code edges()} method
 *  returns its edges.
 *  <p>
 *  This implementation uses <em>Prim's algorithm</em> with an indexed
 *  binary heap.
 *  The constructor takes time proportional to <em>E</em> log <em>V</em>
 *  and extra space (not including the graph) proportional to <em>V</em>,
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the {@code weight()} method takes constant time
 *  and the {@code edges()} method takes time proportional to <em>V</em>.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/43mst">Section 4.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *  For alternate implementations, see {@link LazyPrimMST}, {@link KruskalMST},
 *  and {@link BoruvkaMST}.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class PrimMST implements MST{
  private static final double FLOATING_POINT_EPSILON = 1e-6;
  private IndexMinPQ<Edge> indexMinPQ;
  private boolean[] visisted;
  private double totalWeight;
  private List<Edge> mst;
  //why we need to store all edges incident to the MST?
  //because once the indexed priority deletes an element
  //we cannot retrieve it
  private List<Edge> edgeTo;

  public PrimMST(EdgeWeightedGraph g) {
    assert(g.V()>=0);
    indexMinPQ = new IndexMinPQ<>(g.V());
    visisted = new boolean[g.V()];
    totalWeight = 0;
    mst = new ArrayList<>(g.V());
    edgeTo = new ArrayList<>(g.V());
    for (int i = 0; i < g.V(); i++)
      edgeTo.add(null);
    for (int i = 0; i < g.V(); i++)
      if (!visisted[i])
        prim(g, i);//sometimes there could be a forest of MSTs, not a single MST
  }

  /**
   * starting from v, visit all edges incident
   * to current MST in ascending order of weights
   * dynamically maintain a priority queue of
   * weighted edges for each node adjacent to
   * current MST.
   *
   * @param g
   * @param v
   */
  private void prim(EdgeWeightedGraph g, int v) {
    scan(g, v);
    while(!indexMinPQ.isEmpty()) {
      //we are sure that only node not in current MST
      //will be returned
      int next = indexMinPQ.delMin();
      Edge nextEdge = edgeTo.get(next);
      mst.add(nextEdge);
      totalWeight += nextEdge.weight();
      scan(g, next);
    }
  }

  /**
   * set v to visisted=true
   * update weights of all edges incident to v
   * in the indexed min priority queue
   * @param g
   * @param v
   */
  private void scan(EdgeWeightedGraph g, int v) {
    visisted[v] = true;
    for (Edge e : g.adj(v)) {
      int neighbor = e.other(v);
      if (visisted[neighbor])
        continue;
      if (indexMinPQ.contains(neighbor)) {
        indexMinPQ.decreaseKey(neighbor, e);
        edgeTo.set(neighbor, indexMinPQ.keyOf(neighbor));
      } else {
        indexMinPQ.insert(neighbor, e);
        edgeTo.set(neighbor, e);
      }
    }
  }

  @Override
  public Iterable<Edge> edges() {
    return mst;
  }

  @Override
  public double weight() {
    return totalWeight;
  }
  // check optimality conditions (takes time proportional to E V lg* V)
  private boolean check(EdgeWeightedGraph G) {

    // check weight
    double totalWeight = 0.0;
    for (Edge e : edges()) {
      totalWeight += e.weight();
    }
    if (Math.abs(totalWeight - weight()) > FLOATING_POINT_EPSILON) {
      System.err.printf("Weight of edges does not equal weight(): %f vs. %f\n", totalWeight, weight());
      return false;
    }

    // check that it is acyclic
    UF uf = new UF(G.V());
    for (Edge e : edges()) {
      int v = e.either(), w = e.other(v);
      if (uf.connected(v, w)) {
        System.err.println("Not a forest");
        return false;
      }
      uf.union(v, w);
    }

    // check that it is a spanning forest
    for (Edge e : G.edges()) {
      int v = e.either(), w = e.other(v);
      if (!uf.connected(v, w)) {
        System.err.println("Not a spanning forest");
        return false;
      }
    }

    // check that it is a minimal spanning forest (cut optimality conditions)
    for (Edge e : edges()) {

      // all edges in MST except e
      uf = new UF(G.V());
      for (Edge f : edges()) {
        int x = f.either(), y = f.other(x);
        if (f != e) uf.union(x, y);
      }

      // check that e is min weight edge in crossing cut
      for (Edge f : G.edges()) {
        int x = f.either(), y = f.other(x);
        if (!uf.connected(x, y)) {
          if (f.weight() < e.weight()) {
            System.err.println("Edge " + f + " violates cut optimality conditions");
            return false;
          }
        }
      }

    }

    return true;
  }

  /**
   * Unit tests the {@code PrimMST} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    EdgeWeightedGraph G = new EdgeWeightedGraph(in);
    PrimMST mst = new PrimMST(G);
    for (Edge e : mst.edges()) {
      StdOut.println(e);
    }
    StdOut.printf("%.5f\n", mst.weight());
  }
}
