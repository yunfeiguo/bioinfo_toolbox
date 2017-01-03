package algorithm_practice.leetcode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.*;

public class PutChair {
  public static void main(String[] args) {
    PutChair test = new PutChair();
    System.out.println(test.putChair(new char[][]{new char[]{'E','O','C'},
            new char[]{'C','E','C'},
            new char[]{'C','C','C'}}));
    System.out.println(test.putChair(new char[][]{new char[]{'E'}}));
    System.out.println(test.putChair(new char[][]{new char[]{'C','O','C'},
            new char[]{'C','E','C'},
            new char[]{'C','O','E'}}));
  }
  /**
   convert the matrix to an undirected graph,
   connect adjacent nodes except for obstacles
   do BFS on each equipment and add it up for each
   empty cell, return the cell with minimal distance
   sum
   */
  public List<Integer> putChair(char[][] gym) {
// Write your solution here.
// assume matrix not empty, not null, every E is reachable from every C
    int n = gym.length;
    int m = gym[0].length;
    int[] minDistanceSum = new int[n * m];
    Graph g = new Graph(n * m);
//build graph, skipping obstacles
/*
E O C
|      |
C-E-C
|   |   |
C-C-C


*/
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if (gym[i][j] == 'O') {
          continue;
        }
        //only connect right and down to avoid duplicates
        if (i + 1< n && gym[i + 1][j] != 'O') {
          g.connect(linearize(i, j, m), linearize(i + 1, j, m));
        }
        if (j + 1 < m && gym[i][j + 1] != 'O') {
          g.connect(linearize(i, j, m), linearize(i, j + 1, m));
        }
      }
    }
//calculate min distances by iterating over all equipment
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if(gym[i][j] == 'E') {
          int[] minDistanceToE = g.getDistances(linearize(i, j, m));
          for (int k = 0; k < n * m; k++) {
            minDistanceSum[k] += minDistanceToE[k];
          }
        }
      }
    }
    int indexOfMin = -1; //we can't initialize it to 0, because it might not be a C
    for (int i = 0; i < n * m; i++) {
      List<Integer> coordinates = delinearize(i, m);
      if (gym[coordinates.get(0)][coordinates.get(1)] != 'C') {
        continue;
      }
      if (indexOfMin == -1 || minDistanceSum[i] < minDistanceSum[indexOfMin]) {
        indexOfMin = i;
      }
    }
    return indexOfMin == -1 ? null : delinearize(indexOfMin, m);
  }
  private class Graph {
    /**
     initialize a graph with n points
     */
    private final List<Integer>[] neighbors;
    private final int n;
    public Graph(int n) {
      //assume n > 0
      neighbors = new List[n];
      this.n = n;
      for (int i = 0; i < n; i++) {
        neighbors[i] = new ArrayList<Integer>();
      }
    }
    /**
     use BFS to traverse the entire graph from p
     return shortest distances from p to every node
     */
    public int[] getDistances(int p) {
      int[] distances = new int[n];
      boolean[] visited = new boolean[n];
      Queue<Integer> toVisit = new LinkedList<Integer>();
      int level = 0;
      toVisit.offer(p);
      visited[p] = true;
      while(!toVisit.isEmpty()) {
        int currentLevelSize = toVisit.size();
        for (int i = 0; i < currentLevelSize; i++) {
          int currentNode = toVisit.poll();
          distances[currentNode] = level;
          for (int neighbor : neighbors[currentNode]) {
            if (!visited[neighbor]) {
              toVisit.offer(neighbor);
						/*once it's in the queue, it will be visited for sure
						if mark it visited when setting the level
						some nodes will be repeatedly added
						to queue
						*/
              visited[currentNode] = true;
            }
          }
        }
        level++;
      }
      return distances;
    }
    /**
     build an undirected connection between two nodes
     */
    public void connect(int p, int q) {
      neighbors[p].add(q);
      neighbors[q].add(p);
    }
  }
  /**
   x row index
   y column index
   return linearized index assuming matrix is flattened
   */
  private int linearize(int x, int y, int ncol) {
    return x * ncol + y;
  }
  private List<Integer> delinearize(int i, int ncol) {
    return Arrays.asList(i / ncol, i % ncol);
  }
}
