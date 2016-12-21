package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac Cycle.java
 *  Execution:    java  Cycle filename.txt
 *  Dependencies: Graph.java Stack.java In.java StdOut.java
 *
 *  Identifies a cycle.
 *  Runs in O(E + V) time.
 *
 *  % java Cycle tinyG.txt
 *  3 4 5 3 
 * 
 *  % java Cycle mediumG.txt 
 *  15 0 225 15 
 * 
 *  % java Cycle largeG.txt 
 *  996673 762 840164 4619 785187 194717 996673 
 *
 ******************************************************************************/

import java.util.ArrayList;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>Cycle</tt> class represents a data type for 
 *  determining whether an undirected graph has a cycle.
 *  The <em>hasCycle</em> operation determines whether the graph has
 *  a cycle and, if so, the <em>cycle</em> operation returns one.
 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>hasCycle</em> operation takes constant time;
 *  the <em>cycle</em> operation takes time proportional
 *  to the length of the cycle.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41graph">Section 4.1</a>   
 *  of <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Cycle {
    //use dfs to find cycles in the graph    
    private boolean[] marked;
    private int[] edgeTo;
    private ArrayList<Integer> cycle;
    public Cycle(Graph G) {
        this.cycle = new ArrayList<Integer>();
        if (hasSelfLoop(G)) return;
        if (hasParallelEdges(G)) return;        
        this.marked = new boolean[G.V()];        
        edgeTo = new int[G.V()];
        for (int i = 0; i < G.V(); i++) {
            if (!marked[i]) {
                dfs(G, i, i);
            }
        }        
    }
    private boolean hasSelfLoop(Graph G) {
        for (int v = 0; v < G.V(); v++) {
            for (Integer w : G.adj(v)) {
                if (w == v) {
                    //there is a self-loop
                    cycle.add(v);
                    cycle.add(v);
                    return true;
                }
            }
        }
        return false;
    }
    // does this graph have two parallel edges?
    // side effect: initialize cycle to be two parallel edges (what does it mean??)
    private boolean hasParallelEdges(Graph G) {        
        for (int i = 0; i < G.V(); i++) {
            marked = new boolean[G.V()];            
            for (Integer j : G.adj(i)) {                
                if (marked[j]) {
                    cycle.add(i);
                    cycle.add(j);
                    cycle.add(i);
                    return true;
                } else {
                    marked[j] = true;
                }
            }
        }
        return false;
    }
    public Iterable<Integer> cycle() {
        return cycle;
    }
    private void dfs(Graph G, int v, int vParent) {
        if (hasCycle()) {
            return;
        }
        marked[v] = true;        
        for (Integer i : G.adj(v)) {            
            /*it is possible to go back to immediate ancestor, which has just been visited
             * so we need to know who is the immediate ancestor
             */
            //if (marked[i]) {
            if (i != vParent) {
                if (marked[i]) {
                    for (int x = v; x != i; x = edgeTo[x]) {
                        cycle.add(x);
                    }
                    cycle.add(i);
                    cycle.add(cycle.get(0)); //add the first element to form a cycle
                    return;                
                } else {
                    edgeTo[i] = v;
                    dfs(G, i, v);
                }
            }
            if (hasCycle()) {
                return;
            }
        }
    }
    public boolean hasCycle() {
        return !cycle.isEmpty();
    }
    public static void main(String[] args) {
        In in = new In(args[0]);
        Graph G = new Graph(in);
        Cycle finder = new Cycle(G);
        if (finder.hasCycle()) {
            for (int v : finder.cycle()) {
                StdOut.print(v + " ");
            }
            StdOut.println();
        }
        else {
            StdOut.println("Graph is acyclic");
        }
    }
}

