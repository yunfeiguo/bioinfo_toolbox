package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.SymbolDigraph;
import edu.princeton.cs.algs4.Digraph;

/******************************************************************************
 *  Compilation:  javac Topoological.java
 *  Execution:    java  Topological filename.txt separator
 *  Dependencies: Digraph.java DepthFirstOrder.java DirectedCycle.java
 *                EdgeWeightedDigraph.java EdgeWeightedDirectedCycle.java
 *                SymbolDigraph.java
 *  Data files:   http://algs4.cs.princeton.edu/42digraph/jobs.txt
 *
 *  Compute topological ordering of a DAG or edge-weighted DAG.
 *  Runs in O(E + V) time.
 *
 *  % java Topological jobs.txt "/"
 *  Calculus
 *  Linear Algebra
 *  Introduction to CS
 *  Programming Systems
 *  Algorithms
 *  Theoretical CS
 *  Artificial Intelligence
 *  Machine Learning
 *  Neural Networks
 *  Robotics
 *  Scientific Computing
 *  Computational Biology
 *  Databases
 *
 *
 ******************************************************************************/

/**
 *  The <tt>Topological</tt> class represents a data type for 
 *  determining a topological order of a directed acyclic graph (DAG).
 *  Recall, a digraph has a topological order if and only if it is a DAG.
 *  The <em>hasOrder</em> operation determines whether the digraph has
 *  a topological order, and if so, the <em>order</em> operation
 *  returns one.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>hasOrder</em> and <em>rank</em> operations takes constant time;
 *  the <em>order</em> operation takes time proportional to <em>V</em>.
 *  <p>
 *  See {@link DirectedCycle}, {@link DirectedCycleX}, and
 *  {@link EdgeWeightedDirectedCycle} to compute a
 *  directed cycle if the digraph is not a DAG.
 *  See {@link TopologicalX} for a nonrecursive queue-based algorithm
 *  to compute a topological order of a DAG.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/42digraph">Section 4.2</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Topological {
    private int[] rank;
    private Iterable<Integer> order;
    public Topological(Digraph g) {
        DirectedCycle cycleDetecter = new DirectedCycle(g);
        if (cycleDetecter.hasCycle()) {
            return;
        }
        rank = new int[g.V()];
        DepthFirstOrder dfs = new DepthFirstOrder(g);
        order = dfs.reversePost();
        int i = 0;
        for (int v : order) {
            rank[v] = i++;
        }
    }
    public Topological(EdgeWeightedDigraph g) {
        EdgeWeightedDirectedCycle cycleDetecter = new EdgeWeightedDirectedCycle(g);
        if(cycleDetecter.hasCycle()) {
            throw new IllegalArgumentException("cycle deteced");
        }
        rank = new int[g.V()];
        DepthFirstOrder dfs = new DepthFirstOrder(g);
        order = dfs.reversePost();
        int i = 0;
        for (int v : order) {
            rank[v] = i++;
        }
    }
    /**
     * Returns a topological order if the digraph has a topologial order,
     * and <tt>null</tt> otherwise.
     * @return a topological order of the vertices (as an interable) if the
     *    digraph has a topological order (or equivalently, if the digraph is a DAG),
     *    and <tt>null</tt> otherwise
     */
    public Iterable<Integer> order() {
        if (hasOrder()) {
            return order;
        }
        return null;
    }
    /**
     * Does the digraph have a topological order?
     * @return <tt>true</tt> if the digraph has a topological order (or equivalently,
     *    if the digraph is a DAG), and <tt>false</tt> otherwise
     */
    public boolean hasOrder() {
        return order != null;
    }
    /**
     * The the rank of vertex <tt>v</tt> in the topological order;
     * -1 if the digraph is not a DAG
     * @return the position of vertex <tt>v</tt> in a topological order
     *    of the digraph; -1 if the digraph is not a DAG
     * @throws IndexOutOfBoundsException unless <tt>v</tt> is between 0 and
     *    <em>V</em> &minus; 1
     */
    public int rank(int v) {
        validateVertex(v);
        return rank[v];
    }
    // throw an IndexOutOfBoundsException unless 0 <= v < V
    private void validateVertex(int v) {
        int V = rank.length;
        if (v < 0 || v >= V)
            throw new IndexOutOfBoundsException("vertex " + v + " is not between 0 and " + (V-1));
    }

    /**
     * Unit tests the <tt>Topological</tt> data type.
     */
    public static void main(String[] args) {
        String filename  = args[0];
        String delimiter = args[1];
        SymbolDigraph sg = new SymbolDigraph(filename, delimiter);
        Topological topological = new Topological(sg.G());
        for (int v : topological.order()) {
            StdOut.println(sg.name(v));
        }
    }
}
