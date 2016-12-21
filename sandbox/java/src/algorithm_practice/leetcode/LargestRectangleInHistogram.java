package algorithm_practice.leetcode;

import javax.management.RuntimeErrorException;
import java.util.Arrays;

/**
 * Created by guoy28 on 11/11/16.
 */
public class LargestRectangleInHistogram {
  public int largestRectangleArea(int[] heights) {
    int n = heights.length;
    int maxArea = 0;
    //WeightedQuickUnionFind cluster = new WeightedQuickUnionFind(n);
    WeightedQuickUnionFindWithPathCompression cluster = new WeightedQuickUnionFindWithPathCompression(n);
    Pair[] heightAndIndex = new Pair[n];
    for (int i = 0; i < n; i++) {
      heightAndIndex[i] = new Pair(heights[i], i);
    }
    Arrays.sort(heightAndIndex);
    for (Pair p : heightAndIndex) {
      //we iterate from largest to smallest
      //check left, if saw larger heights, merge with them
      if (p.getIndex() > 0 && heights[p.getIndex() - 1] >= p.getHeight()) {
        cluster.union(p.getIndex() - 1, p.getIndex());
      }
      //check right, ...
      if (p.getIndex() < n - 1 && heights[p.getIndex() + 1] >= p.getHeight()) {
        cluster.union(p.getIndex() + 1, p.getIndex());
      }
      //calculate max area after merging
      maxArea = Math.max(maxArea, p.getHeight() * cluster.size(p.getIndex()));
    }
    return maxArea;
  }
  private class Pair implements Comparable{
    private final int height;
    private final int index;
    public Pair(int h, int i) {
      this.height = h;
      this.index = i;
    }
    public int getHeight() {
      return height;
    }
    public int getIndex() {
      return index;
    }

    @Override
    //we want larger height on top in a minheap
    public int compareTo(Object o) {
      Pair that = (Pair) o;
      return that.height - this.height;
    }
  }

  private class WeightedQuickUnionFindWithPathCompression {
    private int[] roots;
    private int[] sizes;
    private int count;

    public WeightedQuickUnionFindWithPathCompression(int n) {
      roots = new int[n];
      sizes = new int[n];
      count = 0;
      for (int i = 0; i < n; i++) {
        roots[i] = i;
        sizes[i] = 1;
      }
    }
    public int find(int p) {
      while(p != roots[p]) {
        //path compression!
        roots[p] = roots[roots[p]];
        p = roots[p];
      }
      return p;
    }
    //return size of a specific connected component
    public int size(int i) {
      return sizes[find(i)];
    }
    public void union(int p, int q) {
      p = find(p);
      q = find(q);
      if (p == q) {
        return;
      }
      //weighting
      //only expand the larger of two trees
      if(sizes[p] > sizes[q]) {
        roots[q] = p;
        sizes[p] += sizes[q];
      } else {
        roots[p] = q;
        sizes[q] += sizes[p];
      }
      count--;
    }
  }
  private class WeightedQuickUnionFind {
    private int[] roots;
    private int[] sizes;
    private int count;

    public WeightedQuickUnionFind(int n) {
      roots = new int[n];
      sizes = new int[n];
      count = 0;
      for (int i = 0; i < n; i++) {
        roots[i] = i;
        sizes[i] = 1;
      }
    }
    public int find(int p) {
      while(p != roots[p])
        p = roots[p];
      return p;
    }
    //return size of a specific connected component
    public int size(int i) {
      return sizes[find(i)];
    }
    public void union(int p, int q) {
      p = find(p);
      q = find(q);
      if (p == q) {
        return;
      }
      //weighting
      //only expand the larger of two trees
      if(sizes[p] > sizes[q]) {
        roots[q] = p;
        sizes[p] += sizes[q];
      } else {
        roots[p] = q;
        sizes[q] += sizes[p];
      }
      count--;
    }
  }
  public static void main(String[] args) throws Exception {
    LargestRectangleInHistogram test = new LargestRectangleInHistogram();
    if (test.largestRectangleArea(new int[]{2,1,5,6,2,3}) != 10) {
      throw new Exception("wrong answer");
    } else {
      System.out.println("correct!");
    }
  }
}
