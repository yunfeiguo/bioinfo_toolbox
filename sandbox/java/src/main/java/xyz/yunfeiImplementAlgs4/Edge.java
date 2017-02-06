package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/5/17.
 */
/******************************************************************************
 *  Compilation:  javac Edge.java
 *  Execution:    java Edge
 *  Dependencies: StdOut.java
 *
 *  Immutable weighted edge.
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code Edge} class represents a weighted edge in an
 *  {@link EdgeWeightedGraph}. Each edge consists of two integers
 *  (naming the two vertices) and a real-value weight. The data type
 *  provides methods for accessing the two endpoints of the edge and
 *  the weight. The natural order for this data type is by
 *  ascending order of weight.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/43mst">Section 4.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Edge implements Comparable<Edge>{
  private final int s; //one end
  private final int t; //the other end
  private final double w; //weight

  public Edge(int s, int t, double w) {
    if (s < 0 || t < 0)
      throw new IllegalArgumentException("vertex must be non-negative");
    if (Double.isNaN(w))
      throw new IllegalArgumentException("weight is NaN");
    this.s = s;
    this.t = t;
    this.w = w;
  }

  /**
   *
   * @return one end
   */
  public int either() {
    return s;
  }
  public double weight() {
    return w;
  }

  /**
   *
   * @param v one end of the edge
   * @return the other end
   */
  public int other(int v) {
    if (v == s) {
      return t;
    }
    if (v == t) {
      return s;
    }
    throw new IllegalArgumentException(v + ": not an end of this edge");
  }
  @Override
  public int compareTo(Edge o) {
    return Double.compare(w, o.w);
  }
  @Override
  public String toString() {
    return String.format("%d-%d %.5f", s, t, w);
  }
  /**
   * Unit tests the {@code Edge} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    Edge e = new Edge(12, 34, 5.67);
    StdOut.println(e);
  }
}
