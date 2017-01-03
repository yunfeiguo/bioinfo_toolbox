package yunfeiImplementAlgs4;

import java.util.ArrayDeque;
import java.util.Deque;

import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac DepthFirstDirectedPaths.java
 *  Execution:    java DepthFirstDirectedPaths G s
 *  Dependencies: Digraph.java Stack.java
 *
 *  Determine reachability in a digraph from a given vertex using
 *  depth first search.
 *  Runs in O(E + V) time.
 *
 *  % tinyDG.txt 3
 *  3 to 0:  3-5-4-2-0
 *  3 to 1:  3-5-4-2-0-1
 *  3 to 2:  3-5-4-2
 *  3 to 3:  3
 *  3 to 4:  3-5-4
 *  3 to 5:  3-5
 *  3 to 6:  not connected
 *  3 to 7:  not connected
 *  3 to 8:  not connected
 *  3 to 9:  not connected
 *  3 to 10:  not connected
 *  3 to 11:  not connected
 *  3 to 12:  not connected
 *
 ******************************************************************************/

/**
 *  The <tt>DepthFirstDirectedPaths</tt> class represents a data type for finding
 *  directed paths from a source vertex <em>s</em> to every
 *  other vertex in the digraph.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>,
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  It uses extra space (not including the graph) proportional to <em>V</em>.
 *  <p>
 *  For additional documentation,  
 *  see <a href="http://algs4.cs.princeton.edu/42digraph">Section 4.2</a> of  
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne. 
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class DepthFirstDirectedPaths {
    private boolean[] visited;
    private int[] edgeTo; 
    private final int source;
    //private int source; //add final to make it fixed!
    
    public DepthFirstDirectedPaths(Digraph G, int s) {
        visited = new boolean[G.V()];
        edgeTo = new int[G.V()];
        source = s;
        dfs(G, s);
    }
    private void dfs(Digraph G, int v) {
        visited[v] = true;
        for (int w : G.adj(v)) {
            if (!visited[w]) {
                edgeTo[w] = v;
                dfs(G, w);
            }
        }
    }
    public boolean hasPathTo(int v) {
        return visited[v];
    }
    public Iterable<Integer> pathTo(int v) {
        Deque<Integer> path = new ArrayDeque<Integer>();
        if (hasPathTo(v)) {            
            for (int x = v; x != source; x = edgeTo[x]) {
                path.addFirst(x);
            }
            path.addFirst(source);
        }
        return path;
    }
    /**
     * Unit tests the <tt>DepthFirstDirectedPaths</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        // StdOut.println(G);

        int s = Integer.parseInt(args[1]);
        DepthFirstDirectedPaths dfs = new DepthFirstDirectedPaths(G, s);

        for (int v = 0; v < G.V(); v++) {
            if (dfs.hasPathTo(v)) {
                StdOut.printf("%d to %d:  ", s, v);
                for (int x : dfs.pathTo(v)) {
                    if (x == s) StdOut.print(x);
                    else        StdOut.print("-" + x);
                }
                StdOut.println();
            } else {
                StdOut.printf("%d to %d:  not connected\n", s, v);
            }
        }
    }

}
