package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/4/17.
 */

import edu.princeton.cs.algs4.StdOut;

import java.util.Arrays;

/******************************************************************************
 *  Compilation:  javac Biconnected.java
 *  Execution:    java Biconnected V E
 *  Dependencies: Graph.java GraphGenerator.java
 *
 *  Identify articulation points and print them out.
 *  This can be used to decompose a graph into biconnected components.
 *  Runs in O(E + V) time.
 *
 *  an articulation vertex is a vertex, has it being removed, will increase
 *  number of connected components of the graph. a biconnected graph is
 *  a graph that contains no articulation vertex.
 *
 *  http://www.cs.brown.edu/courses/cs016/book/slides/Connectivity2x2.pdf
 *
 ******************************************************************************/
public class Biconnected {
  private int count; //count number of visited vertices in preorder traversal
  private Graph g;
  private int[] preorder; //preorder traversal order of each vertex
  private int[] lowestReachable; //lowest preorder number reachable (in preorder)
  private boolean[] isArticulation;

  public Biconnected(Graph g) {
    assert (g.V() > 0);
    this.g = g;
    preorder = new int[g.V()];
    lowestReachable = new int[g.V()];
    isArticulation = new boolean[g.V()];
    Arrays.fill(preorder, -1);
    for (int i = 0; i < g.V(); i++) {
      if (preorder[i] == -1) {
        dfs(i, i);
      }
    }
  }

  private void dfs(int parent, int node) {
    int childCount = 0;
    preorder[node] = count++;
    lowestReachable[node] = preorder[node];
    for (int neighbor : g.adj(node)) {
      if (preorder[neighbor] == -1) {
        //we only count unvisited children
        //all children in the same cycle will be visited
        //in one recursive call
        childCount++;
        dfs(node, neighbor);
        lowestReachable[node] = Math.min(lowestReachable[node], lowestReachable[neighbor]);
        //in textbook answer, lowestReachable[neighbor] >= preorder[neighbor] was
        //used, I don't think > can happen.
        isArticulation[neighbor] |= lowestReachable[neighbor] == preorder[neighbor];
      } else if (parent != neighbor) {
        lowestReachable[node] = Math.min(lowestReachable[node], preorder[neighbor]);
      }
    }
    //after this node's lowestReachable has been updated
    //we determine if itself is an articulation vertex or not
    isArticulation[node] |= parent == node && childCount > 1;
  }
  public boolean isArticulation(int v) {
    return isArticulation[v];
  }
  public static void main(String[] args) {
    Graph G = new Graph(5);
    G.addEdge(0,1);
    G.addEdge(0,2);
    G.addEdge(0,3);
    G.addEdge(0,4);
    G.addEdge(1,2);
    G.addEdge(3,4);
    StdOut.println(G);

    Biconnected bic = new Biconnected(G);

    // print out articulation points
    StdOut.println();
    StdOut.println("Articulation points");
    StdOut.println("-------------------");
    for (int v = 0; v < G.V(); v++)
      if (bic.isArticulation(v)) StdOut.println(v);

    StdOut.println("2nd graph, suppose no art point");
    Graph G2 = new Graph(5);
    G2.addEdge(0,1);
    G2.addEdge(0,2);
    G2.addEdge(2,1);
    G2.addEdge(3,1);
    G2.addEdge(3,2);
    StdOut.println(G2);

    Biconnected bic2 = new Biconnected(G2);

    // print out articulation points
    StdOut.println();
    StdOut.println("Articulation points");
    StdOut.println("-------------------");
    for (int v = 0; v < G2.V(); v++)
      if (bic2.isArticulation(v)) StdOut.println(v);

    StdOut.println("3rd graph, suppose no art point");
    Graph G3 = new Graph(6);
    G3.addEdge(0,1);
    G3.addEdge(0,4);
    G3.addEdge(0,5);
    G3.addEdge(4,5);
    G3.addEdge(1,2);
    G3.addEdge(3,2);
    G3.addEdge(3,1);
    StdOut.println(G3);

    Biconnected bic3 = new Biconnected(G3);

    // print out articulation points
    StdOut.println();
    StdOut.println("Articulation points");
    StdOut.println("-------------------");
    for (int v = 0; v < G3.V(); v++)
      if (bic3.isArticulation(v)) StdOut.println(v);
  }
}
