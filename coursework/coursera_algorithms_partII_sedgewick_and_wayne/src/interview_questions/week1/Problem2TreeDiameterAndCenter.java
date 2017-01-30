package interview_questions.week1;

import edu.princeton.cs.algs4.BreadthFirstPaths;
import edu.princeton.cs.algs4.Graph;

import java.util.List;

/**
 question:
Diameter and center of a tree. Given a connected graph with no cycles

* Diameter: design a linear-time algorithm to find the longest simple path in
 the graph.
* Center: design a linear-time algorithm to find a vertex such that its maximum
 distance from any other vertex is minimized.

 answer:
 to find diameter, begin with any vertex s, use BFS to find furthest vertex u,
 then from u, use BFS to find furthest vertex v, `$d(u, v) = diameter$`.

 to find center, we look for the middle point on diameter, its maximum distance
 from any other vertex is d/2 or d/2+1.
 */
public class Problem2TreeDiameterAndCenter {
  private Graph g;
  private int diameter;
  private List<Integer> diameterPath;
  private int center;
  public Problem2TreeDiameterAndCenter(Graph g) {
    this.g = g; //should make a copy due to mutability
    this.diameter = 0;
    //assume graph is not empty
    assert(g.V() > 0);
    findDiameterAndCenter();
  }
  private void findDiameterAndCenter() {
    BreadthFirstPaths bfs1 = new BreadthFirstPaths(g, 0);
    int u = 0;
    int maxDistance = 0;
    for (int i = 0; i < g.V(); i++) {
      if (bfs1.distTo(i) > maxDistance) {
        maxDistance = bfs1.distTo(i);
        u = i;
      }
    }
    BreadthFirstPaths bfs2 = new BreadthFirstPaths(g, u);
    int v = 0;
    for (int i = 0; i < g.V(); i++) {
      if (bfs2.distTo(i) > diameter) {
        diameter = bfs2.distTo(i);
        v = i;
      }
    }
    //now locate center
    for (int i : bfs2.pathTo(v)) {
      if (bfs2.distTo(i) >= 0.5*diameter) {
        center = i;
        break;
      }
    }
  }

  public int getDiameter() {
    return diameter;
  }
  public int getCenter() {
    return center;
  }
  public static void main(String[] args) {
    Graph tree = new Graph(10);
    tree.addEdge(0,1);
    tree.addEdge(1,2);
    tree.addEdge(1,3);
    tree.addEdge(3,4);
    tree.addEdge(0,5);
    tree.addEdge(5,6);
    tree.addEdge(5,8);
    tree.addEdge(8,7);
    tree.addEdge(8,9);
    Problem2TreeDiameterAndCenter test = new Problem2TreeDiameterAndCenter(tree);
    System.out.println(test.getDiameter() == 6);
    System.out.println(test.getCenter() == 0);
  }
}

