package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/8/17.
 */
/******************************************************************************
 *  Compilation:  javac FordFulkerson.java
 *  Execution:    java FordFulkerson V E
 *  Dependencies: FlowNetwork.java FlowEdge.java Queue.java
 *  Data files:   http://algs4.cs.princeton.edu/65maxflow/tinyFN.txt
 *
 *  Ford-Fulkerson algorithm for computing a max flow and
 *  a min cut using shortest augmenting path rule.
 *
 java FordFulkerson tinyFN.txt
 6 8
 0:  0->1 0.0/2.0  0->2 0.0/3.0
 1:  1->3 0.0/3.0  1->4 0.0/1.0
 2:  2->3 0.0/1.0  2->4 0.0/1.0
 3:  3->5 0.0/2.0
 4:  4->5 0.0/3.0
 5:

 Max flow from 0 to 5
 0->1 2.0/2.0
 0->2 2.0/3.0
 1->3 1.0/3.0
 1->4 1.0/1.0
 2->3 1.0/1.0
 2->4 1.0/1.0
 3->5 2.0/2.0
 4->5 2.0/3.0
 Min cut: 0 2
 Max flow value = 4.0
 *
 ******************************************************************************/

import com.sun.tools.javac.comp.Flow;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

import java.util.*;

import static experiment.cmdLineParserExamples.Stock.V;
import static yunfeiImplementAlgs4.FlowEdge.FLOATING_POINT_ERROR;

/**
 *  The {@code FordFulkerson} class represents a data type for computing a
 *  <em>maximum st-flow</em> and <em>minimum st-cut</em> in a flow
 *  network.
 *  <p>
 *  This implementation uses the <em>Ford-Fulkerson</em> algorithm with
 *  the <em>shortest augmenting path</em> heuristic.
 *  The constructor takes time proportional to <em>E V</em> (<em>E</em> + <em>V</em>)
 *  in the worst case and extra space (not including the network)
 *  proportional to <em>V</em>, where <em>V</em> is the number of vertices
 *  and <em>E</em> is the number of edges. In practice, the algorithm will
 *  run much faster.
 *  Afterwards, the {@code inCut()} and {@code value()} methods take
 *  constant time.
 *  <p>
 *  If the capacities and initial flow values are all integers, then this
 *  implementation guarantees to compute an integer-valued maximum flow.
 *  If the capacities and floating-point numbers, then floating-point
 *  roundoff error can accumulate.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/64maxflow">Section 6.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class FordFulkerson {
  private final FlowNetwork g;
  private final int V;
  private final int E;
  private final int source;
  private final int sink;
  private double maxflow;
  private List<FlowEdge> edgeTo;
  private boolean[] visited;
  /**
   * find max flow from s to t in flow network g
   * @param g
   * @param s source
   * @param t sink
   */
  public FordFulkerson(FlowNetwork g, int s, int t) {
    this.g = g; //rigorously speaking, we should create a copy of the network
    this.V = g.V();
    this.E = g.E();
    this.source = s;
    this.sink = t;
    this.maxflow = 0;
    this.edgeTo = new ArrayList<>();

    for (int i = 0; i < V; i++)
      edgeTo.add(null);

    while(true) {
      Collection<FlowEdge> path = findAugmentPath();
      if (path == null || path.isEmpty())
        //no more path to augment, quit
        break;
      double bottleneckCapacity = findBottleneckCapacity(path);
      augmentPath(path, bottleneckCapacity);
    }
  }

  /**
   *
   * @return max flow
   */
  public double value() {
    return maxflow;
  }
  /**
   * Returns true if the specified vertex is on the {@code s} side of the mincut.
   *
   * @param  v vertex
   * @return {@code true} if vertex {@code v} is on the {@code s} side of the micut;
   *         {@code false} otherwise
   * @throws IllegalArgumentException unless {@code 0 <= v < V}
   */
  public boolean inCut(int v) {
    validate(v);
    return visited[v];
  }

  /**
   * find bottleneck capacity following path from source to sink
   * @param path
   * @return
   */
  private double findBottleneckCapacity(Collection<FlowEdge> path) {
    double minCapacity = Double.POSITIVE_INFINITY;
    int head = source;
    for (FlowEdge e : path) {
      head = e.other(head);
      minCapacity = Math.min(e.residualCapacityTo(head), minCapacity);
    }
    return minCapacity;
  }

  /**
   * find an augmenting path from source to sink
   * @return
   */
  private Collection<FlowEdge> findAugmentPath() {
    java.util.Queue<Integer> q = new LinkedList<>();
    visited = new boolean[V];
    Deque<FlowEdge> path = new ArrayDeque<>();
    q.offer(source);
    visited[source] = true;
    while(!q.isEmpty()) {
      int next = q.poll();
      for (FlowEdge e : g.adj(next)) {
        int tail = e.other(next);
        if (!visited[tail] && e.residualCapacityTo(tail) > 0) {
          q.offer(tail);
          edgeTo.set(tail, e);
          visited[tail] = true;
          if (tail == sink) {
            for (int i = sink; i != source; i = edgeTo.get(i).other(i)) {
              path.addFirst(edgeTo.get(i));
            }
            return path;
          }
        }
      }
    }
    return path;
  }

  /**
   * apply flow delta onto the path
   * @param path
   */
  private void augmentPath(Collection<FlowEdge> path, double delta) {
    int head = source;
    for (FlowEdge e : path) {
      head = e.other(head);
      e.addResidualFlowTo(head, delta);
    }
    maxflow += delta;
  }
  // throw an IllegalArgumentException if v is outside prescibed range
  private void validate(int v)  {
    if (v < 0 || v >= V)
      throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
  }
  // return excess flow at vertex v
  private double excess(FlowNetwork G, int v) {
    double excess = 0.0;
    for (FlowEdge e : G.adj(v)) {
      if (v == e.from()) excess -= e.flow();
      else               excess += e.flow();
    }
    return excess;
  }

  // return excess flow at vertex v
  private boolean isFeasible(FlowNetwork G, int s, int t) {

    // check that capacity constraints are satisfied
    for (int v = 0; v < G.V(); v++) {
      for (FlowEdge e : G.adj(v)) {
        if (e.flow() < -FLOATING_POINT_ERROR || e.flow() > e.capacity() + FLOATING_POINT_ERROR) {
          System.err.println("Edge does not satisfy capacity constraints: " + e);
          return false;
        }
      }
    }

    // check that net flow into a vertex equals zero, except at source and sink
    if (Math.abs(maxflow + excess(G, s)) > FLOATING_POINT_ERROR) {
      System.err.println("Excess at source = " + excess(G, s));
      System.err.println("Max flow         = " + maxflow);
      return false;
    }
    if (Math.abs(maxflow - excess(G, t)) > FLOATING_POINT_ERROR) {
      System.err.println("Excess at sink   = " + excess(G, t));
      System.err.println("Max flow         = " + maxflow);
      return false;
    }
    for (int v = 0; v < G.V(); v++) {
      if (v == s || v == t) continue;
      else if (Math.abs(excess(G, v)) > FLOATING_POINT_ERROR) {
        System.err.println("Net flow out of " + v + " doesn't equal zero");
        return false;
      }
    }
    return true;
  }
  // check optimality conditions
  private boolean check(FlowNetwork G, int s, int t) {

    // check that flow is feasible
    if (!isFeasible(G, s, t)) {
      System.err.println("Flow is infeasible");
      return false;
    }

    // check that s is on the source side of min cut and that t is not on source side
    if (!inCut(s)) {
      System.err.println("source " + s + " is not on source side of min cut");
      return false;
    }
    if (inCut(t)) {
      System.err.println("sink " + t + " is on source side of min cut");
      return false;
    }

    // check that value of min cut = value of max flow
    double mincutValue = 0.0;
    for (int v = 0; v < G.V(); v++) {
      for (FlowEdge e : G.adj(v)) {
        if ((v == e.from()) && inCut(e.from()) && !inCut(e.to()))
          mincutValue += e.capacity();
      }
    }

    if (Math.abs(mincutValue - value()) > FLOATING_POINT_ERROR) {
      System.err.println("Max flow value = " + value() + ", min cut value = " + mincutValue);
      return false;
    }

    return true;
  }

  /**
   * Unit tests the {@code FordFulkerson} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {

    // create flow network with V vertices and E edges
    In in = new In(args[0]);
    FlowNetwork G = new FlowNetwork(in);
    int s = 0, t = G.V()-1;
    StdOut.println(G);

    // compute maximum flow and minimum cut
    FordFulkerson maxflow = new FordFulkerson(G, s, t);
    StdOut.println("Max flow from " + s + " to " + t);
    for (int v = 0; v < G.V(); v++) {
      for (FlowEdge e : G.adj(v)) {
        if ((v == e.from()) && e.flow() > 0)
          StdOut.println("   " + e);
      }
    }

    // print min-cut
    StdOut.print("Min cut: ");
    for (int v = 0; v < G.V(); v++) {
      if (maxflow.inCut(v)) StdOut.print(v + " ");
    }
    StdOut.println();

    StdOut.println("Max flow value = " +  maxflow.value());
  }
}
