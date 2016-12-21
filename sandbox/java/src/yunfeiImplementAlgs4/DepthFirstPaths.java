package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac DepthFirstPaths.java
 *  Execution:    java DepthFirstPaths G s
 *  Dependencies: Graph.java Stack.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/41graph/tinyCG.txt
 *
 *  Run depth first search on an undirected graph.
 *  Runs in O(E + V) time.
 *
 *  %  java Graph tinyCG.txt
 *  6 8
 *  0: 2 1 5 
 *  1: 0 2 
 *  2: 0 1 3 4 
 *  3: 5 4 2 
 *  4: 3 2 
 *  5: 3 0 
 *
 *  % java DepthFirstPaths tinyCG.txt 0
 *  0 to 0:  0
 *  0 to 1:  0-2-1
 *  0 to 2:  0-2
 *  0 to 3:  0-2-3
 *  0 to 4:  0-2-3-4
 *  0 to 5:  0-2-3-5
 *
 ******************************************************************************/

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>DepthFirstPaths</tt> class represents a data type for finding
 *  paths from a source vertex <em>s</em> to every other vertex
 *  in an undirected graph.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>,
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  It uses extra space (not including the graph) proportional to <em>V</em>.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41graph">Section 4.1</a>   
 *  of <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class DepthFirstPaths {
    private int[] edgeTo; //edgeTo[v] = last edge on s-v path
    private boolean[] marked; //marked[v] = is there a path from s to v?
    private final int s; //source vertex
    
    /**
     * computes a path between <tt>s</tt> and every other vertex in graph G
     * @param G
     * @param s
     */
    public DepthFirstPaths(Graph G, int s) {
        edgeTo = new int[G.V()];
        marked = new boolean[G.V()];
        this.s = s;
        dfs(G, s);
    }
    private void dfs(Graph G, int v) {
        marked[v] = true;
        for (int w : G.adj(v)) {
            if (!marked[w]) {
                edgeTo[w] = v;
                dfs(G, w);
            }
        }
    }
    /**
     * is there a path between source s and vertex v
     * @param v
     * @return
     */
    public boolean hasPathTo(int v) {
        return marked[v];
    }
    /**
     * 
     * @param v
     * @return the sequence of vertices on a path between the source vertex 
     * and the vertex, as an Iterable
     */
    public Iterable<Integer> pathTo(int v) {
        Deque<Integer> path = new ArrayDeque<Integer>();
        if (!hasPathTo(v)) {
            return path;
        }
        for (int w = v; w != s; w = edgeTo[w]) {
            path.addFirst(w);
        }
        path.addFirst(s);
        return path;
    }
    /**
     * Unit tests the <tt>DepthFirstPaths</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Graph G = new Graph(in);
        int s = Integer.parseInt(args[1]);
        DepthFirstPaths dfs = new DepthFirstPaths(G, s);

        for (int v = 0; v < G.V(); v++) {
            if (dfs.hasPathTo(v)) {
                StdOut.printf("%d to %d:  ", s, v);
                for (int x : dfs.pathTo(v)) {
                    if (x == s) StdOut.print(x);
                    else        StdOut.print("-" + x);
                }
                StdOut.println();
            }

            else {
                StdOut.printf("%d to %d:  not connected\n", s, v);
            }

        }
    }
}
