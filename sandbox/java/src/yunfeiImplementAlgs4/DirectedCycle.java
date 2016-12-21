package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac DirectedCycle.java
 *  Execution:    java DirectedCycle input.txt
 *  Dependencies: Digraph.java Stack.java StdOut.java In.java
 *  Data files:   http://algs4.cs.princeton.edu/42digraph/tinyDG.txt
 *                http://algs4.cs.princeton.edu/42digraph/tinyDAG.txt
 *
 *  Finds a directed cycle in a digraph.
 *  Runs in O(E + V) time.
 *
 *  % java DirectedCycle tinyDG.txt 
 *  Directed cycle: 3 5 4 3 
 *
 *  %  java DirectedCycle tinyDAG.txt 
 *  No directed cycle
 *
 ******************************************************************************/

/**
 *  The <tt>DirectedCycle</tt> class represents a data type for 
 *  determining whether a digraph has a directed cycle.
 *  The <em>hasCycle</em> operation determines whether the digraph has
 *  a directed cycle and, and of so, the <em>cycle</em> operation
 *  returns one.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>hasCycle</em> operation takes constant time;
 *  the <em>cycle</em> operation takes time proportional
 *  to the length of the cycle.
 *  <p>
 *  See {@link Topological} to compute a topological order if the
 *  digraph is acyclic.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/42digraph">Section 4.2</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
import java.util.ArrayList;
import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

public class DirectedCycle {
    private boolean[] visited;
    private boolean[] onStack; //is an element on the path being visited
    private ArrayList<Integer> cycle;
    private int[] edgeTo;
    public DirectedCycle(Digraph g) {
        visited = new boolean[g.V()];
        edgeTo = new int[g.V()];
        onStack = new boolean[g.V()];
        for (int i = 0; i < g.V(); i++) {
            if (!visited[i]) {
                dfs(g, i);    
            }
        }
    }
    private void dfs(Digraph g, int v) {
        /*
        if (hasCycle())
            return; //think about whether two hasCycle() are necessary
        */
        visited[v] = true;
        onStack[v] = true;
        for (int i : g.adj(v)) {
            if (hasCycle()) return;
            if (onStack[i]) {
                //there is a cycle
                cycle = new ArrayList<Integer>();
                cycle.add(i);
                for (int p = v; p != i; p = edgeTo[p])
                    cycle.add(p);
                cycle.add(i);
                return;
            } else if (!visited[i]){
                edgeTo[i] = v;
                dfs(g, i);
                /*
                if (hasCycle())                 
                    return;
                //think about whether two hasCycle() are necessary
                 */                
            }
        }
        onStack[v] = false;
    }
    boolean hasCycle() {
        return cycle != null;
    }
    public Iterable<Integer> cycle() {
        if (hasCycle()) {
            return cycle;
        }
        return null;
    }
    /**
     * Unit tests the <tt>DirectedCycle</tt> data type.
     */
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);

        DirectedCycle finder = new DirectedCycle(G);
        if (finder.hasCycle()) {
            StdOut.print("Directed cycle: ");
            for (int v : finder.cycle()) {
                StdOut.print(v + " ");
            }
            StdOut.println();
        }

        else {
            StdOut.println("No directed cycle");
        }
        StdOut.println();
    }
}
