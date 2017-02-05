package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/4/17.
 */

import edu.princeton.cs.algs4.StdOut;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

/******************************************************************************
 *  Compilation:  javac AllPaths.java
 *  Execution:    java AllPaths
 *  Depedencies:  Graph.java
 *
 *  Enumerate all simple paths (of length >= 1) in a graph between s and t.
 *  This implementation uses depth-first search and backtracking.
 *
 *  Warning: there can be exponentially many simple paths in a graph,
 *           so no algorithm is suitable for large graphs.
 *
 *  7 vertices, 9 edges
 *  0: 2 1
 *  1: 5 0
 *  2: 5 3 0
 *  3: 6 4 2
 *  4: 6 5 3
 *  5: 4 1 2
 *  6: 4 3
 *
 *
 *  all simple paths between 0 and 6:
 *  0-2-5-4-6
 *  0-2-5-4-3-6
 *  0-2-3-6
 *  0-2-3-4-6
 *  0-1-5-4-6
 *  0-1-5-4-3-6
 *  0-1-5-2-3-6
 *  0-1-5-2-3-4-6
 *  # paths = 8
 *
 *  all simple paths between 1 and 5:
 *  1-5
 *  1-0-2-5
 *  1-0-2-3-6-4-5
 *  1-0-2-3-4-5
 *  # paths = 4
 *
 ******************************************************************************/
public class AllPaths {
  private List<List> paths;
  private Deque<Integer> verticesOnPath;
  private Graph g;
  private boolean[] marked;
  /**
   * find all simple paths from s to t in graph g
   * @param g
   * @param s
   * @param t
   */
  public AllPaths(Graph g, int s, int t) {
    this.g = g;
    paths = new ArrayList<>();
    verticesOnPath = new ArrayDeque<>();
    marked = new boolean[g.V()];
    assert(g.degree(s) != 0);
    assert(g.degree(t) != 0);
    assert(s != t);
    verticesOnPath.addFirst(s);
    marked[s] = true;
    dfs(s, t);
  }

  /**
   * dfs traverse all simple paths that end at t
   * @param s
   * @param t
   */
  private void dfs(int s, int t) {
    if (s == t) {
      paths.add(new ArrayList(verticesOnPath));
      return;
    }
    for (int i : g.adj(s)) {
      if (marked[i])
        continue;
      verticesOnPath.addFirst(i);
      marked[i] = true;
      dfs(i, t);
      verticesOnPath.removeFirst();
      marked[i] = false;
    }
  }
  public int numberOfPaths() {
    return paths.size();
  }
  // test client
  public static void main(String[] args) {
    Graph G = new Graph(7);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(2, 5);
    G.addEdge(1, 5);
    G.addEdge(5, 4);
    G.addEdge(3, 6);
    G.addEdge(4, 6);
    StdOut.println(G);

    StdOut.println();
    StdOut.println("all simple paths between 0 and 6:");
    AllPaths allpaths1 = new AllPaths(G, 0, 6);
    StdOut.println("# paths = " + allpaths1.numberOfPaths());

    StdOut.println();
    StdOut.println("all simple paths between 1 and 5:");
    AllPaths allpaths2 = new AllPaths(G, 1, 5);
    StdOut.println("# paths = " + allpaths2.numberOfPaths());
  }
}
