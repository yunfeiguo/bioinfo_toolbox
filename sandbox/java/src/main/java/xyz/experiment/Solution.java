package experiment;

/**
 * Created by guoy28 on 11/5/16.
 */

import java.util.List;
import java.util.*;

public class Solution {
  public List<Integer> closest(int[] a, int[] b, int[] c, int k) {
    // Write your solution here.
    //assume a,b,c not null, not empty, and there are enough elements for finding kth point
    Set<String> seen = new HashSet<String>(); //assume no overflow
    Queue<Point> q = new PriorityQueue<Point>(k, new PointComparator());
    Point first = new Point(0, 0, 0);
    seen.add(first.toString());
    q.offer(first);
    for (int i = 1; i < k; i++) {
      Point closest = q.poll();
      Point nextA = closest.copy();
      nextA.increaseA();
      Point nextB = closest.copy();
      nextB.increaseB();
      Point nextC = closest.copy();
      nextC.increaseC();
      //we cannot simply skip points with previously seen distances because
      //we may miss some candidates
      if (seen.add(nextA.toString())) {
        q.offer(nextA);
      }
      if (seen.add(nextB.toString())) {
        q.offer(nextB);
      }
      if (seen.add(nextC.toString())) {
        q.offer(nextC);
      }
    }
    return q.peek().toList(a, b, c);
  }
  //private class PointComparator<Point> implements Comparator<Point> {
  private class PointComparator implements Comparator<Point> {
    private int[] a;
    private int[] b;
    private int[] c;
    public PointComparator() {}
    public PointComparator(int x) {}
    public boolean equals(Object obj) {
      return false; //no implementation, should throw exception
    }
    public int compare(Point x, Point y) {
      return 0;
    }
  }
  private class Point {
    private int first; //indices to underlying int array
    private int second;
    private int third;
    public Point(int one, int two, int three) {
      first = one;
      second = two;
      third = three;
    }
    public Point copy() {
      return new Point(first, second, third);
    }
    public void increaseA() {
      this.first++;
    }
    public void increaseB() {
      this.second++;
    }
    public void increaseC() {
      this.third++;
    }
    public List<Integer> toList(int[] a, int[] b, int[] c) {
      List<Integer> result = new ArrayList<Integer>();
      result.add(a[first]);
      result.add(b[second]);
      result.add(c[third]);
      return result;
    }
    /**
     calculate Euclidean distance given array indices and arrays
     return 0 if any index out of bound
     */
    public long distance(int[] a, int[] b, int[] c) {
      long distance = 0;
      if (first >= a.length || second >= b.length || third >= c.length) {
        return 0;
      }
      distance += (long) a[first] * a[first];
      distance += (long) b[second] * b[second];
      distance += (long) c[third] * c[third];
      return distance;
    }
    @Override
    public String toString() {
      //why space? prevent 11 1 <-> 1 11
      return first + " " + second + " " + third;
    }
  }

  public static void main(String[] args) {
    Solution test = new Solution();
    test.closest(new int[]{1, 2}, new int[]{1,2}, new int[]{1,2}, 1);
  }
}

