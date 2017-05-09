package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/24/17.
 */
/******************************************************************************
 *  Compilation:  javac DirectedEdge.java
 *  Execution:    java DirectedEdge
 *  Dependencies: StdOut.java
 *
 *  Immutable weighted directed edge.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code DirectedEdge} class represents a weighted edge in an
 *  {@link EdgeWeightedDigraph}. Each edge consists of two integers
 *  (naming the two vertices) and a real-value weight. The data type
 *  provides methods for accessing the two endpoints of the directed edge and
 *  the weight.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/44sp">Section 4.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class DirectedEdge {
  private final int from;
  private final int to;
  private final double weight;

  public DirectedEdge(int from, int to, double weight) {
    if (from < 0 || to < 0) throw new IllegalArgumentException("cannot be negative");
    if (Double.isNaN(weight)) throw new IllegalArgumentException("cannot be NaN");
    this.from = from;
    this.to = to;
    this.weight = weight;
  }
  public int from() {
    return from;
  }
  public int to() {
    return to;
  }
  public double weight() {
    return weight;
  }
  @Override
  public String toString() {
    return from + "->" + to + " " + String.format("%5.2f", weight);
  }
  /**
   * Unit tests the {@code DirectedEdge} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    DirectedEdge e = new DirectedEdge(12, 34, 5.67);
    StdOut.println(e);
  }
}
