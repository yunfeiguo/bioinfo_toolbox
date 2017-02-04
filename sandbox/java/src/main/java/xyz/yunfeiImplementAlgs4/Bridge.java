package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/2/17.
 */

import edu.princeton.cs.algs4.StdOut;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/******************************************************************************
 *  Compilation:  javac Bridge.java
 *  Execution:    java  Bridge V E
 *  Dependencies: Graph.java GraphGenerator.java
 *
 *  Identifies bridge edges and prints them out. This decomposes
 *  a directed graph into two-edge connected components.
 *  Runs in O(E + V) time.
 *
 *  Key quantity:  low[v] = minimum DFS preorder number of v
 *  and the set of vertices w for which there is a back edge (x, w)
 *  with x a descendant of v and w an ancestor of v.
 *
 *  Note: code assumes no parallel edges, e.g., two parallel edges
 *  would be (incorrectly) identified as bridges.
 *
 ******************************************************************************/
public class Bridge {
  private int count;
  private List<int[]> bridges;
  private int[] preorder; //preorder sequence
  private int[] lowestReachable; //lowest node reachable in preorder traversal
  private Graph g;
  public Bridge(Graph g) {
    assert(g.V()> 0);
    this.g = g;
    bridges = new ArrayList<>();
    preorder = new int[g.V()];
    lowestReachable = new int[g.V()];
    count = 0;
    Arrays.fill(preorder, -1);
    for (int i = 0; i < g.V(); i++) {
      if (preorder[i] == -1) {//unvisited
        dfs(i, i);
      }
    }
  }

  /**
   * DFS traversal of child
   * here we pass parent such that we do not
   * look back
   *
   * during traversal, mark preorder sequence
   * and identify bridge by comparing lowestReachable
   * and preorder
   *
   * @param parent
   * @param child
   */
  private void dfs(int parent, int child) {
    preorder[child] = count++;
    lowestReachable[child] = preorder[child];
    for (int neighbor : g.adj(child)) {
      if (preorder[neighbor] == -1) {
        dfs(child, neighbor);
        //updating what child can reach traversing DOWN DFS tree
        lowestReachable[child] = Math.min(lowestReachable[child], lowestReachable[neighbor]);
        if (lowestReachable[neighbor] == preorder[neighbor]) {
          bridges.add(new int[]{child, neighbor});
        }
      } else if (parent != neighbor) {
        //we do not look back
        //if we have a neighbor which has been visited before, then we must
        //be in a cycle
        lowestReachable[child] = Math.min(lowestReachable[child], preorder[neighbor]);
      }
    }
  }

  /**
   *
   * @return how many connected components are bridged?
   */
  public int components() {
    return bridges.size() + 1;
  }

  /**
   *
   * @return iterable of pairs of nodes at each bridge
   */
  public Iterable<int[]> getBridges() {
    return new ArrayList(bridges);
  }
  // test client
  public static void main(String[] args) {
    Graph G = new Graph(6);
    G.addEdge(0,1);
    G.addEdge(0,2);
    G.addEdge(2,1);
    G.addEdge(3,1);
    G.addEdge(3,4);
    G.addEdge(4,5);
    G.addEdge(3,5);
    StdOut.println(G);

    Bridge bridge = new Bridge(G);
    StdOut.println("Edge connected components = " + bridge.components());
    for (int[] b : bridge.getBridges()) {
      if (b.length == 2)
        StdOut.println("bridge " + b[0] + " " + b[1]);
    }
    StdOut.println("another example");
    Graph G2 = new Graph(4);
    G2.addEdge(0,1);
    G2.addEdge(0,2);
    G2.addEdge(2,1);
    G2.addEdge(3,1);
    G2.addEdge(3,2);
    StdOut.println(G2);

    Bridge bridge2 = new Bridge(G2);
    StdOut.println("Edge connected components = " + bridge2.components());
    for (int[] b : bridge2.getBridges()) {
      if (b.length == 2)
        StdOut.println("bridge " + b[0] + " " + b[1]);
    }
  }
}
