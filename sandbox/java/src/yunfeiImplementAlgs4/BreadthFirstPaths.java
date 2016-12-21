package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac BreadthFirstPaths.java
 *  Execution:    java BreadthFirstPaths G s
 *  Dependencies: Graph.java Queue.java Stack.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/41graph/tinyCG.txt
 *
 *  Run breadth first search on an undirected graph.
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
 *  %  java BreadthFirstPaths tinyCG.txt 0
 *  0 to 0 (0):  0
 *  0 to 1 (1):  0-1
 *  0 to 2 (1):  0-2
 *  0 to 3 (2):  0-2-3
 *  0 to 4 (2):  0-2-4
 *  0 to 5 (1):  0-5
 *
 *  %  java BreadthFirstPaths largeG.txt 0
 *  0 to 0 (0):  0
 *  0 to 1 (418):  0-932942-474885-82707-879889-971961-...
 *  0 to 2 (323):  0-460790-53370-594358-780059-287921-...
 *  0 to 3 (168):  0-713461-75230-953125-568284-350405-...
 *  0 to 4 (144):  0-460790-53370-310931-440226-380102-...
 *  0 to 5 (566):  0-932942-474885-82707-879889-971961-...
 *  0 to 6 (349):  0-932942-474885-82707-879889-971961-...
 *
 ******************************************************************************/

import java.util.ArrayDeque;
import java.util.Deque;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>BreadthFirstPaths</tt> class represents a data type for finding
 *  shortest paths (number of edges) from a source vertex <em>s</em>
 *  (or a set of source vertices)
 *  to every other vertex in an undirected graph.
 *  <p>
 *  This implementation uses breadth-first search.
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
public class BreadthFirstPaths {
    private static final int INFINITY = Integer.MAX_VALUE;
    private boolean[] marked;
    private int[] edgeTo;
    private int[] distTo;
    /**
     * Computes the shortest path between the source vertex <tt>s</tt>
     * and every other vertex in the graph <tt>G</tt>.
     * @param G the graph
     * @param s the source vertex
     */
    public BreadthFirstPaths(Graph G, int s) {
        marked = new boolean[G.V()];
        edgeTo = new int[G.V()];
        distTo = new int[G.V()];
        bfs(G, s);

        assert check(G, s);
    }
    /**
     * Computes the shortest path between any one of the source vertices in <tt>sources</tt>
     * and every other vertex in graph <tt>G</tt>.
     * @param G the graph
     * @param sources the source vertices
     */
    public BreadthFirstPaths(Graph G, Iterable<Integer> sources) {
        marked = new boolean[G.V()];
        distTo = new int[G.V()];
        edgeTo = new int[G.V()];
        bfs(G, sources);
    }

    private void bfs(Graph G, int s) {
        Queue<Integer> toVisit = new Queue<Integer>();
        for (int v = 0; v < G.V(); v++) {
            distTo[v] = INFINITY;
        }
        toVisit.enqueue(s);
        marked[s] = true;
        distTo[s] = 0;
        while (!toVisit.isEmpty()) {
            int current = toVisit.dequeue();
            for (int i : G.adj(current)) {
                if (!marked[i]) {
                    marked[i] = true; //we must mark it before putting the vertex into queue, to avoid repeated visits.
                    edgeTo[i] = current;
                    distTo[i] = distTo[current] + 1;
                    toVisit.enqueue(i);
                }
            }
        }
    }
    private void bfs(Graph G, Iterable<Integer> sources) {
        Queue<Integer> toVisit = new Queue<Integer>();
        for (int v = 0; v < G.V(); v++) {
            distTo[v] = INFINITY;
        }
        for (int s : sources) {                    
            toVisit.enqueue(s);
            marked[s] = true;
            distTo[s] = 0;
        }
        while (!toVisit.isEmpty()) {
            int current = toVisit.dequeue();
            for (int i : G.adj(current)) {
                if (!marked[i]) {
                    marked[i] = true;
                    edgeTo[i] = current;
                    distTo[i] = distTo[current] + 1;
                    toVisit.enqueue(i);
                }
            }
        }
    }
    public boolean hasPathTo(int v) {
        return marked[v];
    }
    public int distTo(int v) {
        return distTo[v];
    }
    /**
     * Returns a shortest path between the source vertex <tt>s</tt> (or sources)
     * and <tt>v</tt>, or <tt>null</tt> if no such path.
     * @param v the vertex
     * @return the sequence of vertices on a shortest path, as an Iterable
     */
    public Iterable<Integer> pathTo(int v) {
        Deque<Integer> path = new ArrayDeque<Integer>();
        if (!hasPathTo(v)) {
            return path;
        }
        int i;
        for (i = v; distTo[i] != 0; i = edgeTo[i] ) {
            path.addFirst(i);
        }
        path.addFirst(i);        
        return path;
    }
    // check optimality conditions for single source
    private boolean check(Graph G, int s) {

        // check that the distance of s = 0
        if (distTo[s] != 0) {
            StdOut.println("distance of source " + s + " to itself = " + distTo[s]);
            return false;
        }

        // check that for each edge v-w dist[w] <= dist[v] + 1
        // provided v is reachable from s
        for (int v = 0; v < G.V(); v++) {
            for (int w : G.adj(v)) {
                if (hasPathTo(v) != hasPathTo(w)) {
                    StdOut.println("edge " + v + "-" + w);
                    StdOut.println("hasPathTo(" + v + ") = " + hasPathTo(v));
                    StdOut.println("hasPathTo(" + w + ") = " + hasPathTo(w));
                    return false;
                }
                if (hasPathTo(v) && (distTo[w] > distTo[v] + 1)) {
                    StdOut.println("edge " + v + "-" + w);
                    StdOut.println("distTo[" + v + "] = " + distTo[v]);
                    StdOut.println("distTo[" + w + "] = " + distTo[w]);
                    return false;
                }
            }
        }

        // check that v = edgeTo[w] satisfies distTo[w] = distTo[v] + 1
        // provided v is reachable from s
        for (int w = 0; w < G.V(); w++) {
            if (!hasPathTo(w) || w == s) continue;
            int v = edgeTo[w];
            if (distTo[w] != distTo[v] + 1) {
                StdOut.println("shortest path edge " + v + "-" + w);
                StdOut.println("distTo[" + v + "] = " + distTo[v]);
                StdOut.println("distTo[" + w + "] = " + distTo[w]);
                return false;
            }
        }

        return true;
    }

    /**
     * Unit tests the <tt>BreadthFirstPaths</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Graph G = new Graph(in);
        // StdOut.println(G);

        int s = Integer.parseInt(args[1]);
        int count = 0;
        BreadthFirstPaths bfs = new BreadthFirstPaths(G, s);

        for (int v = 0; v < G.V(); v++) {
            if (bfs.hasPathTo(v)) {
                StdOut.printf("%d to %d (%d):  ", s, v, bfs.distTo(v));
                for (int x : bfs.pathTo(v)) {
                    if (x == s) StdOut.print(x);
                    else        StdOut.print("-" + x);
                }
                StdOut.println();
                count++;
                if (count > 5) {
                    break;
                }
            }

            else {
                StdOut.printf("%d to %d (-):  not connected\n", s, v);
            }           

        }
    }
}
