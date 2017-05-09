package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac DepthFirstOrder.java
 *  Execution:    java DepthFirstOrder filename.txt
 *  Dependencies: Digraph.java Queue.java Stack.java StdOut.java
 *                EdgeWeightedDigraph.java DirectedEdge.java
 *  Data files:   http://algs4.cs.princeton.edu/42digraph/tinyDAG.txt
 *                http://algs4.cs.princeton.edu/42digraph/tinyDG.txt
 *
 *  Compute preorder and postorder for a digraph or edge-weighted digraph.
 *  Runs in O(E + V) time.
 *
 *  % java DepthFirstOrder tinyDAG.txt
 *     v  pre post
 *  --------------
 *     0    0    8
 *     1    3    2
 *     2    9   10
 *     3   10    9
 *     4    2    0
 *     5    1    1
 *     6    4    7
 *     7   11   11
 *     8   12   12
 *     9    5    6
 *    10    8    5
 *    11    6    4
 *    12    7    3
 *  Preorder:  0 5 4 1 6 9 11 12 10 2 3 7 8 
 *  Postorder: 4 5 1 12 11 10 9 6 0 3 2 7 8 
 *  Reverse postorder: 8 7 2 3 0 6 9 10 11 12 1 5 4 
 *
 ******************************************************************************/

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Stack;

import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>DepthFirstOrder</tt> class represents a data type for 
 *  determining depth-first search ordering of the vertices in a digraph
 *  or edge-weighted digraph, including preorder, postorder, and reverse postorder.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>preorder</em>, <em>postorder</em>, and <em>reverse postorder</em>
 *  operation takes take time proportional to <em>V</em>.
 *  <p>
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/42digraph">Section 4.2</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class DepthFirstOrder {
    private ArrayList<Integer> preorder;
    private ArrayList<Integer> postorder;
    private Deque<Integer> reversePostOrder;
    private int preRank;
    private int postRank;
    private int[] pre;
    private int[] post;
    private boolean[] visited;
    
    public DepthFirstOrder(Digraph g) {
        //assume it's DAG
        preRank = 0;
        postRank = 0;
        preorder = new ArrayList<Integer>();
        postorder = new ArrayList<Integer>();
        reversePostOrder = new ArrayDeque<Integer>();
        pre = new int[g.V()];
        post = new int[g.V()];
        visited = new boolean[g.V()];
        for (int i = 0; i < g.V(); i++) {
            if (!visited[i]) {
                dfs(g, i);
            }
        }
    }
    public DepthFirstOrder(EdgeWeightedDigraph g) {
        //assume it's DAG
        preRank = 0;
        postRank = 0;
        preorder = new ArrayList<Integer>();
        postorder = new ArrayList<Integer>();
        reversePostOrder = new ArrayDeque<Integer>();
        pre = new int[g.V()];
        post = new int[g.V()];
        visited = new boolean[g.V()];
        for (int i = 0; i < g.V(); i++) {
            if (!visited[i]) {
                dfs(g, i);
            }
        }
    }
    private void dfs(EdgeWeightedDigraph g, int v) {
        visited[v] = true;
        preorder.add(v);
        pre[v] = preRank++;
        for (DirectedEdge e: g.adj(v)) {
            int to = e.to();
            if (!visited[to]) {
                dfs(g, to);
            }
        }
        postorder.add(v);
        reversePostOrder.addFirst(v);
        post[v] = postRank++;
    }
    private void dfs(Digraph g, int v) {
        visited[v] = true;
        preorder.add(v);
        pre[v] = preRank++;
        for (int i : g.adj(v)) {
            if (!visited[i]) {
                dfs(g, i);
            }
        }
        postorder.add(v);
        reversePostOrder.addFirst(v);
        post[v] = postRank++;
    }
    
    /**
     * Returns the preorder number of vertex <tt>v</tt>.
     * @param v the vertex
     * @return the preorder number of vertex <tt>v</tt>
     */
    public int pre(int v) {
        return pre[v];
    }

    /**
     * Returns the postorder number of vertex <tt>v</tt>.
     * @param v the vertex
     * @return the postorder number of vertex <tt>v</tt>
     */
    public int post(int v) {
        return post[v];
    }

    /**
     * Returns the vertices in postorder.
     * @return the vertices in postorder, as an iterable of vertices
     */
    public Iterable<Integer> post() {
        return postorder;
    }

    /**
     * Returns the vertices in preorder.
     * @return the vertices in preorder, as an iterable of vertices
     */
    public Iterable<Integer> pre() {
        return preorder;
    }

    /**
     * Returns the vertices in reverse postorder.
     * @return the vertices in reverse postorder, as an iterable of vertices
     */
    public Iterable<Integer> reversePost() {
        //alternative method
        /*ArrayList<Integer> reverse = new ArrayList<Integer>(postorder);
        Collections.reverse(reverse);
        return reverse;*/
        return reversePostOrder;
    }
    // throw an IllegalArgumentException unless {@code 0 <= v < V}
    private void validateVertex(int v) {
        int V = visited.length;
        if (v < 0 || v >= V)
            throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
    }
    // check that pre() and post() are consistent with pre(v) and post(v)
    private boolean check() {

        // check that post(v) is consistent with post()
        int r = 0;
        for (int v : post()) {
            if (post(v) != r) {
                StdOut.println("post(v) and post() inconsistent");
                return false;
            }
            r++;
        }

        // check that pre(v) is consistent with pre()
        r = 0;
        for (int v : pre()) {
            if (pre(v) != r) {
                StdOut.println("pre(v) and pre() inconsistent");
                return false;
            }
            r++;
        }


        return true;
    }

    /**
     * Unit tests the <tt>DepthFirstOrder</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);

        DepthFirstOrder dfs = new DepthFirstOrder(G);
        StdOut.println("   v  pre post");
        StdOut.println("--------------");
        for (int v = 0; v < G.V(); v++) {
            StdOut.printf("%4d %4d %4d\n", v, dfs.pre(v), dfs.post(v));
        }

        StdOut.print("Preorder:  ");
        for (int v : dfs.pre()) {
            StdOut.print(v + " ");
        }
        StdOut.println();

        StdOut.print("Postorder: ");
        for (int v : dfs.post()) {
            StdOut.print(v + " ");
        }
        StdOut.println();

        StdOut.print("Reverse postorder: ");
        for (int v : dfs.reversePost()) {
            StdOut.print(v + " ");
        }
        StdOut.println();


    }

}
