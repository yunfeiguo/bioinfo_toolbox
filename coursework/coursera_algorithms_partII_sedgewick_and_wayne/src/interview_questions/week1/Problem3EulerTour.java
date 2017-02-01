package interview_questions.week1;

import edu.princeton.cs.algs4.DepthFirstSearch;
import edu.princeton.cs.algs4.Graph;

import javax.management.RuntimeErrorException;
import java.util.*;

/**question:
 *
# Euler cycle. An Euler cycle in a graph is a cycle (not necessarily simple) that
 uses every edge in the graph exactly once.

* Show that a connected graph has an Euler cycle if and only if every vertex has
 * even degree.
* Design a linear-time algorithm to determine whether a graph has an Euler cycle,
 * and if so, find one.

 answer:
Suppose we have a graph with just 2 vertices A and B, connected by 2n edges,
 n = 1,2,3,... then we can always go A->B->A to traverse 2 edges. We can repeat this
 n times to traverse all edges exactly once.

Now suppose we have k vertices in a connected graph with even degree and there is
 an Euler tour in it. Now we add a vertex V and 2n (n=1,2,...,k) edges, each edge connecting V and
 a vertex in the graph. It is easy to see that we can still construct an Euler tour
 by adding k loops to the previous Euler tour.

 if there is a vertex with odd degrees, it is obvious that a cycle cannot be found
 using all the edges from that vertex.

 therefore, an Euler cycle exists on a connected graph if and only if all vertices
 have even degrees.

 but how to quickly find an Euler cycle?
 use DFS but when there are multiple candidates for next stop,
 choose only the one with more than 1 degree left, unless there
 is no other choice (in this case, we are simply coming back
 to the starting point). time complexity: O(E), space complexity
 O(V + E).
 */
public class Problem3EulerTour {
  private Graph g;
  private HashGraph hg;
  private List<Integer> path;
  private boolean isEulerTourExists;
  private int[] edgeTo;
  private List<Integer> tour;

  public boolean isEulerTourExists() {
    //determine if graph is connected
    //determine if degree of every vertex is even
    DepthFirstSearch dfs = new DepthFirstSearch(g, 0);
    if (dfs.count() != g.V())
      return false;
    for (int i = 0; i < g.V(); i++) {
      if (g.degree(i) % 2 != 0)
        return false;
    }
    return true;
  }

  public Problem3EulerTour(Graph g) {
    this.g = g;
    path = new LinkedList<Integer>();
    hg = new HashGraph(g.V());
    //assume graph not empty
    for (int i = 0; i < g.V(); i++) {
      for (int j : g.adj(i)) {
        if (i <= j)
          //make sure no duplicates
          hg.addEdge(i, j);
      }
    }
    isEulerTourExists = isEulerTourExists();
    if (!isEulerTourExists)
      return;
    //use modified dfs to find an Euler tour
    dfs(0);
  }

  /**
   * this graph is designed to change dynamically,
   * i.e., edges can be removed.
   * accordingly a map data structure is used
   * to record counts of edges.
   */
  private class HashGraph extends Graph {
    //record who are neighbors and how many edges are between them
    private Map[] neighbors;
    public HashGraph(int V) {
      super(V);
      neighbors = new Map[V];
      for (int i = 0; i < V; i++) {
        neighbors[i] = new HashMap<Integer, Integer>();
      }
    }
    @Override
    public int degree(int i) {
      return neighbors[i].keySet().size();
    }

    /**
     * add edge between i and j
     * increase count for both i and j
     * @param i
     * @param j
     */
    @Override
    public void addEdge(int i, int j) {
      if (neighbors[i].containsKey(j)) {
        neighbors[i].put(j, ((Map<Integer, Integer>) neighbors[i]).get(j) + 1);
        neighbors[j].put(i, ((Map<Integer, Integer>) neighbors[j]).get(i) + 1);
      } else {
        neighbors[i].put(j, 1);
        neighbors[j].put(i, 1);
      }
    }
    public void removeEdge(int i, int j) {
      if (neighbors[i].containsKey(j)) {
        neighbors[i].put(j, ((Map<Integer, Integer>) neighbors[i]).get(j) - 1);
        neighbors[j].put(i, ((Map<Integer, Integer>) neighbors[i]).get(j));
        if (((Map<Integer, Integer>) neighbors[i]).get(j) <= 0) {
          neighbors[i].remove(j);
          neighbors[j].remove(i);
        }
      }
    }
    @Override
    public Iterable<Integer> adj(int i) {
      return neighbors[i].keySet();
    }
  }

  public Collection<Integer> getPath() {
    return path;
  }

  /**
   * time complexity should be O(E) as each edge is visited once
   * @param source
   */
  private void dfs(int source) {
    //for each vertex, use a HashMap to
    //record which edge has been not visited/visited.
    path.add(source);
    for (int current = source; hg.degree(current) > 0;) {
      int next = -1; //first choice for next step
      int backupNext = -1; //second choice for next step
      for (int neighbor : hg.adj(current)) {
        if (hg.degree(neighbor) > 1) {
          next = neighbor;
          break;
        } else if (hg.degree(neighbor) > 0) {
          //this statement will only be executed when we come back
          // to the starting point, closing the loop
          backupNext = neighbor;
        }
      }
      if (next == -1 && backupNext == -1) {
        throw new RuntimeException("no next step available at node ");
      }
      if (next == -1)
        next = backupNext;
      path.add(next);
      hg.removeEdge(current, next);
      current = next;
    }
  }

  public static void main(String[] args) {
    Graph g = new Graph(6);
    g.addEdge(0,1);
    g.addEdge(0,3);
    g.addEdge(1,3);
    g.addEdge(4,3);
    g.addEdge(1,2);
    g.addEdge(2,3);
    g.addEdge(1,5);
    g.addEdge(5,4);
    Problem3EulerTour test = new Problem3EulerTour(g);
    System.out.println(test.getPath());
  }
}
