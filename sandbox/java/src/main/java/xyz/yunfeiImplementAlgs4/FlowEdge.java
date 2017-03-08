package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/7/17.
 */
/******************************************************************************
 *  Compilation:  javac FlowEdge.java
 *  Execution:    java FlowEdge
 *  Dependencies: StdOut.java
 *
 *  Capacitated edge with a flow in a flow network.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code FlowEdge} class represents a capacitated edge with a
 * flow in a {@link FlowNetwork}. Each edge consists of two integers
 *  (naming the two vertices), a real-valued capacity, and a real-valued
 *  flow. The data type provides methods for accessing the two endpoints
 *  of the directed edge and the weight. It also provides methods for
 *  changing the amount of flow on the edge and determining the residual
 *  capacity of the edge.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/64maxflow">Section 6.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class FlowEdge {
  private final int from;
  private final int to;
  private final double capacity;
  private double flow;
  /**
   * Initializes an edge from vertex {@code v} to vertex {@code w} with
   * the given {@code capacity} and zero flow.
   * @param v the tail vertex
   * @param w the head vertex
   * @param capacity the capacity of the edge
   * @throws IllegalArgumentException if either {@code v} or {@code w}
   *    is a negative integer
   * @throws IllegalArgumentException if {@code capacity < 0.0}
   */
  public FlowEdge(int v, int w, double capacity) {
    if (v < 0) throw new IllegalArgumentException("vertex index must be a non-negative integer");
    if (w < 0) throw new IllegalArgumentException("vertex index must be a non-negative integer");
    if (!(capacity >= 0.0)) throw new IllegalArgumentException("Edge capacity must be non-negative");
    from = v;
    to = w;
    this.capacity = capacity;
  }
  /**
   * Initializes an edge from vertex {@code v} to vertex {@code w} with
   * the given {@code capacity} and {@code flow}.
   * @param v the tail vertex
   * @param w the head vertex
   * @param capacity the capacity of the edge
   * @param flow the flow on the edge
   * @throws IllegalArgumentException if either {@code v} or {@code w}
   *    is a negative integer
   * @throws IllegalArgumentException if {@code capacity} is negative
   * @throws IllegalArgumentException unless {@code flow} is between
   *    {@code 0.0} and {@code capacity}.
   */
  public FlowEdge(int v, int w, double capacity, double flow) {
    this(v, w, capacity);
    if (!(flow <= capacity)) throw new IllegalArgumentException("flow exceeds capacity");
    if (!(flow >= 0.0))      throw new IllegalArgumentException("flow must be non-negative");
    this.flow = flow;
  }
  public FlowEdge(FlowEdge e) {
    this(e.from, e.to, e.capacity, e.flow);
  }
  public int from() {
    return from;
  }
  public int to() {
    return to;
  }
  /**
   * Returns the endpoint of the edge that is different from the given vertex
   * (unless the edge represents a self-loop in which case it returns the same vertex).
   * @param v one endpoint of the edge
   * @return the endpoint of the edge that is different from the given vertex
   *   (unless the edge represents a self-loop in which case it returns the same vertex)
   * @throws IllegalArgumentException if {@code vertex} is not one of the endpoints
   *   of the edge
   */
  public int other(int v) {
    if (v == from)
      return to;
    else if (v == to)
      return from;
    else
      throw new IllegalArgumentException("not an end");
  }
  /**
   * Returns the residual capacity of the edge in the direction
   *  to the given {@code vertex}.
   * @param v one endpoint of the edge
   * @return the residual capacity of the edge in the direction to the given vertex
   *   If {@code vertex} is the tail vertex, the residual capacity equals
   *   {@code capacity() - flow()}; if {@code vertex} is the head vertex, the
   *   residual capacity equals {@code flow()}.
   * @throws IllegalArgumentException if {@code vertex} is not one of the endpoints of the edge
   */
  public double residualCapacityTo(int v) {
    if (v == from)
      return flow;
    else if (v == to)
      return capacity - flow;
    else
      throw new IllegalArgumentException("not an end");
  }
  /**
   * Increases the flow on the edge in the direction to the given vertex.
   *   If {@code vertex} is the tail vertex, this increases the flow on the edge by {@code delta};
   *   if {@code vertex} is the head vertex, this decreases the flow on the edge by {@code delta}.
   * @param v one endpoint of the edge
   * @param delta amount by which to increase flow
   * @throws IllegalArgumentException if {@code vertex} is not one of the endpoints
   *   of the edge
   * @throws IllegalArgumentException if {@code delta} makes the flow on
   *   on the edge either negative or larger than its capacity
   * @throws IllegalArgumentException if {@code delta} is {@code NaN} or negative
   */
  public void addResidualFlowTo(int v, double delta) {
    if (Double.isNaN(delta) || delta < 0)
      throw new IllegalArgumentException("delta is NaN or negative");
    if  (v != from && v != to)
      throw new IllegalArgumentException("not an end");
    if (v == from) {
      flow -= delta;
    } else if (v == to) {
      flow += delta;
    }
    //here we cannot assume delta is positive, therefore, check flow at the end
    if (flow < 0)
      throw new IllegalArgumentException("flow is negative");
    if (flow > capacity) {
      throw new IllegalArgumentException("flow cannot be larger than capacity");
    }
  }

  @Override
  public String toString() {
    return from + "->" + to + " " + flow + "/" + capacity;
  }
  /**
   * Unit tests the {@code FlowEdge} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    FlowEdge e = new FlowEdge(12, 23, 4.56);
    StdOut.println(e);
  }
}
