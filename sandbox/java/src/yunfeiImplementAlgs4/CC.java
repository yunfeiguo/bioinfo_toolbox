package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac CC.java
 *  Execution:    java CC filename.txt
 *  Dependencies: Graph.java StdOut.java Queue.java
 *  Data files:   http://algs4.cs.princeton.edu/41graph/tinyG.txt
 *
 *  Compute connected components using depth first search.
 *  Runs in O(E + V) time.
 *
 *  % java CC tinyG.txt
 *  3 components
 *  0 1 2 3 4 5 6
 *  7 8 
 *  9 10 11 12
 *
 *  % java CC mediumG.txt 
 *  1 components
 *  0 1 2 3 4 5 6 7 8 9 10 ...
 *
 *  % java -Xss50m CC largeG.txt 
 *  1 components
 *  0 1 2 3 4 5 6 7 8 9 10 ...
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>CC</tt> class represents a data type for 
 *  determining the connected components in an undirected graph.
 *  The <em>id</em> operation determines in which connected component
 *  a given vertex lies; the <em>connected</em> operation
 *  determines whether two vertices are in the same connected component;
 *  the <em>count</em> operation determines the number of connected
 *  components; and the <em>size</em> operation determines the number
 *  of vertices in the connect component containing a given vertex.

 *  The <em>component identifier</em> of a connected component is one of the
 *  vertices in the connected component: two vertices have the same component
 *  identifier if and only if they are in the same connected component.

 *  <p>
 *  This implementation uses depth-first search.
 *  The constructor takes time proportional to <em>V</em> + <em>E</em>
 *  (in the worst case),
 *  where <em>V</em> is the number of vertices and <em>E</em> is the number of edges.
 *  Afterwards, the <em>id</em>, <em>count</em>, <em>connected</em>,
 *  and <em>size</em> operations take constant time.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41graph">Section 4.1</a>   
 *  of <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class CC {
    private int[] id; //record connected component assignment for each vertex
    private boolean[] visited;
    private int[] size; //record # of vertices in each CC
    private int ccCount;
    public CC(Graph G) {
        this.id = new int[G.V()];
        this.visited = new boolean[G.V()];
        this.size = new int[G.V()];
        this.ccCount = 0;
        for (int i = 0; i < G.V(); i++) {
            if (!this.visited[i]) {
                size[ccCount] = 0;
                dfs(G, i);
                ccCount++;
            }            
        }
    }
    private void dfs(Graph G, int v) {        
        visited[v] = true;
        id[v] = ccCount;
        size[ccCount]++;
        for (int u : G.adj(v)) {
            if (!visited[u]) {
                dfs(G, u);
            }
        }        
    }
    public int id(int v) {
        return id[v];
    }
    public int size(int v) {
        return size[id[v]];
    }
    public int count() {
        return ccCount;
    }
    public boolean connected(int v, int w) {
        return id[v] == id[w];
    }
    /**
     * @deprecated
     * @param v
     * @param w
     * @return
     */
    public boolean areConnected(int v, int w) {
        return false;
    }
    public static void main(String[] args) {
        In in = new In(args[0]);
        Graph G = new Graph(in);
        CC cc = new CC(G);

        // number of connected components
        int M = cc.count();
        StdOut.println(M + " components");

        // compute list of vertices in each connected component
        @SuppressWarnings("unchecked")
        Queue<Integer>[] components = (Queue<Integer>[]) new Queue[M];
        for (int i = 0; i < M; i++) {
            components[i] = new Queue<Integer>();
        }
        for (int v = 0; v < G.V(); v++) {
            components[cc.id(v)].enqueue(v);
        }

        // print results
        for (int i = 0; i < M; i++) {
            for (int v : components[i]) {
                StdOut.print(v + " ");
            }
            StdOut.println();
        }
    }

}
