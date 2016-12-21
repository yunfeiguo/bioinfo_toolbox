package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac Bipartite.java
 *  Execution:    java  Bipartite V E F
 *  Dependencies: Graph.java 
 *
 *  Given a graph, find either (i) a bipartition or (ii) an odd-length cycle.
 *  Runs in O(E + V) time.
 *
 *
 ******************************************************************************/

import java.util.ArrayList;

import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.StdRandom;

/**
 *  The <tt>Bipartite</tt> class represents a data type for 
 *  determining whether an undirected graph is bipartite or whether
 *  it has an odd-length cycle.
 *  The <em>isBipartite</em> operation determines whether the graph is
 *  bipartite. If so, the <em>color</em> operation determines a
 *  bipartition; if not, the <em>oddCycle</em> operation determines a
 *  cycle with an odd number of edges.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>isBipartite</em> and <em>color</em> operations
 *  take constant time; the <em>oddCycle</em> operation takes time proportional
 *  to the length of the cycle.
 *  See {@link BipartiteX} for a nonrecursive version that uses breadth-first
 *  search.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41graph">Section 4.1</a>   
 *  of <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Bipartite {
    private boolean isBipartite;
    private ArrayList<Integer> cycle;
    private boolean[] color;
    private boolean[] visited;
    private int[] edgeTo;
    private int V;

    /**
     * Determines whether an undirected graph is bipartite and finds either a
     * bipartition or an odd-length cycle.
     *
     * @param  G the graph
     */
    public Bipartite(Graph G) {
        if(hasSelfLoop(G)) return;
        this.isBipartite = true;
        this.cycle = null;
        this.color = new boolean[G.V()];
        this.visited = new boolean[G.V()];
        this.edgeTo = new int[G.V()];
        this.V = G.V();

        for (int i = 0; i < G.V(); i++) {
            if (!this.visited[i]) {
                color[i] = true;
                dfs(G, i, i);
            }
        }
        assert check(G);
    }
    private boolean check(Graph G) {
        // graph is bipartite
        if (isBipartite) {
            for (int v = 0; v < G.V(); v++) {
                for (int w : G.adj(v)) {
                    if (color[v] == color[w]) {
                        System.err.printf("edge %d-%d with %d and %d in same side of bipartition\n", v, w, v, w);
                        return false;
                    }
                }
            }
        }

        // graph has an odd-length cycle
        else {
            // verify cycle
            int first = -1, last = -1;
            for (int v : oddCycle()) {
                if (first == -1) first = v;
                last = v;
            }
            if (first != last) {
                System.err.printf("cycle begins with %d and ends with %d\n", first, last);
                return false;
            }
        }

        return true;
    }
    private boolean hasSelfLoop(Graph G) {
        for (int v = 0; v < G.V(); v++) {
            for (Integer w : G.adj(v)) {
                if (w == v) {
                    //there is a self-loop
                    cycle = new ArrayList<Integer>();
                    cycle.add(v);
                    cycle.add(v);
                    isBipartite = false;
                    return true;
                }
            }
        }
        return false;
    }
    /**
     * use dfs to traverse the graph, ignore any 
     * immediate ancestor; if encounter a visited node
     * then we can check its color, if different color
     * then bipartite still holds, otherwise, there is 
     * an odd-length cycle
     * @param G
     * @param v
     * @param vParent immediate parent of v
     */
    private void dfs(Graph G, int v, int vParent) {
        if (!isBipartite()) {
            return;
        }
        visited[v] = true;
        for (int i : G.adj(v)) {
            if (i != vParent) {
                if (visited[i] && color(v) == color(i)) {
                    //we got a cycle, and two adjacent nodes with same color
                    //this must be an odd cycle                    
                    isBipartite = false;
                    cycle = new ArrayList<Integer>();
                    for (int x = v; x != i; x = edgeTo[x]) {
                        cycle.add(x);
                    }
                    cycle.add(i);
                    cycle.add(v);
                    return;
                } else if (!visited[i]) {
                    edgeTo[i] = v;
                    color[i] = !color[v];
                    dfs(G, i, v);
                    if (!isBipartite()) {
                        return;
                    }
                }
            }
        }
    }
    /**
     * Returns true if the graph is bipartite.
     *
     * @return <tt>true</tt> if the graph is bipartite; <tt>false</tt> otherwise
     */
    public boolean isBipartite() {
        return isBipartite;
    }
    /**
     * Returns the side of the bipartite that vertex <tt>v</tt> is on.
     *
     * @param  v the vertex
     * @return the side of the bipartition that vertex <tt>v</tt> is on; two vertices
     *         are in the same side of the bipartition if and only if they have the
     *         same color
     * @throws IllegalArgumentException unless <tt>0 &le; v &lt; V</tt> 
     * @throws UnsupportedOperationException if this method is called when the graph
     *         is not bipartite
     */
    public boolean color(int v) {
        if (v < 0 || v > V) {
            throw new IllegalArgumentException();
        }
        if (!isBipartite()) {
            throw new UnsupportedOperationException();
        }
        return color[v];
    }
    /**
     * Returns an odd-length cycle if the graph is not bipartite, and
     * <tt>null</tt> otherwise.
     *
     * @return an odd-length cycle if the graph is not bipartite
     *         (and hence has an odd-length cycle), and <tt>null</tt>
     *         otherwise
     */
    public Iterable<Integer> oddCycle() {
        return cycle;
    }
    /**
     * Unit tests the <tt>Bipartite</tt> data type.
     */
    public static void main(String[] args) {
        int V1 = Integer.parseInt(args[0]);
        int V2 = Integer.parseInt(args[1]);
        int E  = Integer.parseInt(args[2]);
        int F  = Integer.parseInt(args[3]);

        // create random bipartite graph with V1 vertices on left side,
        // V2 vertices on right side, and E edges; then add F random edges
        Graph G = GraphGenerator.bipartite(V1, V2, E);
        for (int i = 0; i < F; i++) {
            int v = StdRandom.uniform(V1 + V2);
            int w = StdRandom.uniform(V1 + V2);
            G.addEdge(v, w);
        }

        StdOut.println(G);


        Bipartite b = new Bipartite(G);
        if (b.isBipartite()) {
            StdOut.println("Graph is bipartite");
            for (int v = 0; v < G.V(); v++) {
                StdOut.println(v + ": " + b.color(v));
            }
        }
        else {
            StdOut.print("Graph has an odd-length cycle: ");
            for (int x : b.oddCycle()) {
                StdOut.print(x + " ");
            }
            StdOut.println();
        }
    }
}
