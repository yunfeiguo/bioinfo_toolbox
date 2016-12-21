package yunfeiImplementAlgs4;

import java.util.Queue;

import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

import java.util.ArrayDeque;
/******************************************************************************
 *  Compilation:  javac BreadthFirstDirectedPaths.java
 *  Execution:    java BreadthFirstDirectedPaths V E
 *  Dependencies: Digraph.java Queue.java Stack.java
 *
 *  Run breadth first search on a digraph.
 *  Runs in O(E + V) time.
 *
 *  % java BreadthFirstDirectedPaths tinyDG.txt 3
 *  3 to 0 (2):  3->2->0
 *  3 to 1 (3):  3->2->0->1
 *  3 to 2 (1):  3->2
 *  3 to 3 (0):  3
 *  3 to 4 (2):  3->5->4
 *  3 to 5 (1):  3->5
 *  3 to 6 (-):  not connected
 *  3 to 7 (-):  not connected
 *  3 to 8 (-):  not connected
 *  3 to 9 (-):  not connected
 *  3 to 10 (-):  not connected
 *  3 to 11 (-):  not connected
 *  3 to 12 (-):  not connected
 *
 ******************************************************************************/

/**
 *  The <tt>BreadthDirectedFirstPaths</tt> class represents a data type for finding
 *  shortest paths (number of edges) from a source vertex <em>s</em>
 *  (or set of source vertices) to every other vertex in the digraph.
 *  <p>
 *  This implementation uses breadth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>,
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  It uses extra space (not including the digraph) proportional to <em>V</em>.
 *  <p>
 *  For additional documentation, 
 *  see <a href="http://algs4.cs.princeton.edu/42digraph">Section 4.2</a> of 
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class BreadthFirstDirectedPaths {
    private static final int INFINITY = Integer.MAX_VALUE;
    private int[] edgeTo;    
    private boolean[] visited;
    private int[] distTo;

    public BreadthFirstDirectedPaths(Digraph g, int s) {
        edgeTo = new int[g.V()];
        visited = new boolean[g.V()];
        distTo = new int[g.V()];
        for (int i = 0; i < g.V(); i++) {
            distTo[i] = INFINITY;
        }
        bfs(g, s);              
    }
    public BreadthFirstDirectedPaths(Digraph g, Iterable<Integer> sources) {
        edgeTo = new int[g.V()];
        visited = new boolean[g.V()];
        distTo = new int[g.V()];
        for (int i = 0; i < g.V(); i++) {
            distTo[i] = INFINITY;
        }
        bfs(g, sources);
    }
    private void bfs(Digraph g, int v) {
        Queue<Integer> q = new ArrayDeque<Integer>();
        distTo[v] = 0;
        visited[v] = true;
        q.offer(v);
        while (!q.isEmpty()) {
            int current = q.poll();
            for (int i : g.adj(current)) {
                if(!visited[i]) {
                    edgeTo[i] = current;
                    visited[i] = true;
                    distTo[i] = distTo[current] + 1;
                    q.offer(i);
                }
            }
        }        
    }
    private void bfs(Digraph g, Iterable<Integer> sources) {
        Queue<Integer> q = new ArrayDeque<Integer>();
        for (int v : sources) {
            distTo[v] = 0;
            visited[v] = true;
            q.offer(v);
        }
        while (!q.isEmpty()) {
            int current = q.poll();
            for (int i : g.adj(current)) {
                if (!visited[i]) {
                    edgeTo[i] = current;
                    visited[i] = true;
                    distTo[i] = distTo[current] + 1;
                    q.offer(i);
                }
            }
        }        
    }    
    public boolean hasPathTo(int v) {
        return visited[v];
    }
    public int distTo(int v) {
        return distTo[v];
    }
    public Iterable<Integer> pathTo(int v) {
        ArrayDeque<Integer> path = new ArrayDeque<Integer>();
        if (hasPathTo(v)) {
            int x = v;
            for (;distTo[x] != 0; x = edgeTo[x]) {
                path.addFirst(x);
            }
            path.addFirst(x);
        }
        return path;
    }
    /**
     * Unit tests the <tt>BreadthFirstDirectedPaths</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        // StdOut.println(G);

        int s = Integer.parseInt(args[1]);
        BreadthFirstDirectedPaths bfs = new BreadthFirstDirectedPaths(G, s);

        for (int v = 0; v < G.V(); v++) {
            if (bfs.hasPathTo(v)) {
                StdOut.printf("%d to %d (%d):  ", s, v, bfs.distTo(v));
                for (int x : bfs.pathTo(v)) {
                    if (x == s) StdOut.print(x);
                    else        StdOut.print("->" + x);
                }
                StdOut.println();
            }

            else {
                StdOut.printf("%d to %d (-):  not connected\n", s, v);
            }

        }
    }

}
