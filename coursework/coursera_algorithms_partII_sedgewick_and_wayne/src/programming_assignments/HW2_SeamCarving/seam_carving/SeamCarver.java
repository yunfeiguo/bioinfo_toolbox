package programming_assignments.HW2_SeamCarving.seam_carving;

import edu.princeton.cs.algs4.Picture;

import java.awt.Color;
import java.util.*;

/**
 * Created by guoy28 on 3/4/17.
 */
public class SeamCarver {
  private static int BORDER_ENERGY = 1000;
  private Picture pic;
  private double[][] energy;
  // create a seam carver object based on the given picture
  public SeamCarver(Picture picture) {
    updatePicture(picture);
  }

  /**
   * replace old picture with new one
   * update energy matrix
   * @param p
   */
  private void updatePicture(Picture p) {
    pic = new Picture(p);
    energy = new double[width()][height()];
    for (int x = 0; x < width(); x++) {
      for (int y = 0; y < height(); y++) {
        energy[x][y] = energy(x, y);
      }
    }
  }
  // current picture
  public Picture picture() {
    return new Picture(pic);
  }
  // width of current picture
  public int width() {
    return pic.width();
  }
  // height of current picture
  public int height() {
    return pic.height();
  }
  // energy of pixel at column x and row y
  public double energy(int x, int y) {
    if (x < 0 || x >= width() || y < 0 || y >= height())
      throw new IndexOutOfBoundsException("coordinates out of bound");
    if (x == 0 || x == width() - 1 || y == 0 || y == height() - 1)
      return BORDER_ENERGY;
    Color left = pic.get(x - 1, y);
    Color right = pic.get(x + 1, y);
    Color top = pic.get(x, y - 1);
    Color bottom = pic.get(x, y + 1);
    double rX = left.getRed() - right.getRed();
    double gX = left.getGreen() - right.getGreen();
    double bX = left.getBlue() - right.getBlue();
    double rY = top.getRed() - bottom.getRed();
    double gY = top.getGreen() - bottom.getGreen();
    double bY = top.getBlue() - bottom.getBlue();
    double deltaX2 = rX * rX + gX * gX + bX * bX;
    double deltaY2 = rY * rY + gY * gY + bY * bY;
    return Math.sqrt(deltaX2 + deltaY2);
  }

  /**
   *
   * sequence of indices for horizontal seam
   * seam has minimum energy among all possible
   * seams
   *
   * **********************
   * **********************
   * --------------------->
   * **********************
   * **********************
   *
    */
  public int[] findHorizontalSeam() {
    NeighborFinder neighborFinder = new HorizontalNeighborFinder(height(), width());
    double minEnergy = Double.POSITIVE_INFINITY;
    int minEnergyRow = -1;
    List<Point> firstColumn = new ArrayList<>(width());
    for (int y = 0; y < height(); y++) {
      firstColumn.add(new Point(0, y));
    }
    MatrixDigraphDAGSP sp = new MatrixDigraphDAGSP(energy, neighborFinder, firstColumn);
    for (int y = 0; y < height(); y++) {
      if (sp.distTo(new Point(width() - 1, y)) < minEnergy) {
        minEnergy = sp.distTo(new Point(width() - 1, y));
        minEnergyRow = y;
      }
    }
    int[] seam = new int[width()];
    int i = 0;
    for (Point p : sp.pathTo(new Point(width() - 1, minEnergyRow))) {
      seam[i++] = p.getY();
    }
    return seam;
  }

  /** sequence of indices for vertical seam
   * seam has minimum energy among all possible
   * seams
   *
   * ********|*************
   * ********|*************
   * ********|*************
   * ********|*************
   * ********\*************
   */
  public int[] findVerticalSeam() {
    NeighborFinder neighborFinder = new VerticalNeighborFinder(height(), width());
    double minEnergy = Double.POSITIVE_INFINITY;
    int minEnergyColumn = -1;
    List<Point> firstRow = new ArrayList<>(width());
    for (int x = 0; x < width(); x++) {
      firstRow.add(new Point(x, 0));
    }
    MatrixDigraphDAGSP sp = new MatrixDigraphDAGSP(energy, neighborFinder, firstRow);
    for (int x = 0; x < width(); x++) {
      if (sp.distTo(new Point(x, height() - 1)) < minEnergy) {
        minEnergy = sp.distTo(new Point(x, height() - 1));
        minEnergyColumn = x;
      }
    }
    int[] seam = new int[height()];
    int i = 0;
    for (Point p : sp.pathTo(new Point(minEnergyColumn, height() - 1))) {
      seam[i++] = p.getX();
    }
    return seam;
  }

  /**
   * remove horizontal seam from current picture
   * @param seam
   * @throws
   */
  public void removeHorizontalSeam(int[] seam) {
    if (seam == null)
      throw new NullPointerException("nothing to remove");
    if (seam.length != width())
      throw new IllegalArgumentException("unequal to width");
    if (height() <= 1)
      throw new IllegalArgumentException("image too small");
    for (int i : seam) {
      if (i < 0 || i >= height())
        throw new IllegalArgumentException("seam index out of bound");
    }
    for (int i = 1; i < seam.length - 1; i++) {
      if (Math.abs(seam[i] - seam[i - 1]) > 1 || Math.abs(seam[i + 1] - seam[i]) > 1)
        throw new IllegalArgumentException("neighbors not differ <= 1");
    }
    if (width() > 1 && (Math.abs(seam[1] - seam[0]) > 1 || Math.abs(seam[width() - 1] - seam[width() - 2]) > 1))
      throw new IllegalArgumentException("neighbors not differ <= 1");
    Picture newPic = new Picture(width(), height() - 1);
    for (int x = 0; x < width(); x++) {
      for (int y = 0; y < height(); y++) {
        if (seam[x] > y)
          newPic.set(x, y, pic.get(x, y));
        else if (seam[x] < y)
          newPic.set(x, y - 1, pic.get(x, y));
      }
    }
    updatePicture(newPic);
  }

  /** remove vertical seam from current picture
   *
   * @param seam
   */
  public void removeVerticalSeam(int[] seam) {
    if (seam == null)
      throw new NullPointerException("nothing to remove");
    if (seam.length != height())
      throw new IllegalArgumentException("unequal to height");
    if (width() <= 1)
      throw new IllegalArgumentException("image too small");
    for (int i : seam) {
      if (i < 0 || i >= width())
        throw new IllegalArgumentException("seam index out of bound");
    }
    for (int i = 1; i < seam.length - 1; i++) {
      if (Math.abs(seam[i] - seam[i - 1]) > 1 || Math.abs(seam[i + 1] - seam[i]) > 1)
        throw new IllegalArgumentException("neighbors not differ <= 1");
    }
    if (height() > 1 && (Math.abs(seam[1] - seam[0]) > 1 || Math.abs(seam[height() - 1] - seam[height() - 2]) > 1))
      throw new IllegalArgumentException("neighbors not differ <= 1");

    Picture newPic = new Picture(width() - 1, height());
    for (int x = 0; x < width(); x++) {
      for (int y = 0; y < height(); y++) {
        if (seam[y] > x)
          newPic.set(x, y, pic.get(x, y));
        else if (seam[y] < x)
          newPic.set(x - 1, y, pic.get(x, y));
      }
    }
    updatePicture(newPic);
  }

  public static void main(String[] args) {
  }
  /**
   * matrix-based directed graph
   * edges defined dynamically
   * find shortest path (SP) from sources
   * to all possible destinations
   */
  private class MatrixDigraphDAGSP {
    private double[][] matrix;
    //right now we use wxh matrix to store distances
    //memory usage can be reduced to linear
    private double[][] distTo;
    private Point[][] closestParentOf;

    /**
     * constructor will find shortest paths from sources
     * to all possible destinations
     *
     * since the edges are dynamically defined by the
     * neighborfinder, we do not really know whether
     * the graph is a DAG or not. To make this class more
     * widely applicable, we need to check its prerequisites.
     * For now, ignore it.
     *
     * @param matrix
     * @param neighborFinder
     * @param sources
     */
    public MatrixDigraphDAGSP(double[][] matrix, NeighborFinder neighborFinder, Collection<Point> sources) {
      if (matrix == null || matrix[0] == null)
        throw new NullPointerException("matrix null");
      this.matrix = matrix;
      int w = matrix.length;
      int h = matrix[0].length;
      if (w == 0)
        throw new IllegalArgumentException("matrix empty");
      distTo = new double[w][h];
      closestParentOf = new Point[w][h];
      //initialization
      for (int x = 0; x < w; x++) {
        Arrays.fill(distTo[x], Double.POSITIVE_INFINITY);
        Arrays.fill(closestParentOf[x], null);
      }
      for (Point p : sources) {
        distTo[p.getX()][p.getY()] = matrix[p.getX()][p.getY()];
      }
      //assume matrix is DAG
      Queue<Point> q = new ArrayDeque<>(sources);
      while(!q.isEmpty()) {
        Set<Point> nextSetInTopologicalOrder = new HashSet<>();
        while (!q.isEmpty()) {
          Point next = q.poll();
          for (Point p : neighborFinder.adj(next)) {
            relax(next, p);
            nextSetInTopologicalOrder.add(p);
          }
        }
        //we can do this because we know the graph is being
        //visited in topological order
        for(Point p : nextSetInTopologicalOrder) {
          q.offer(p);
        }
      }
    }

    /**
     * return
     * @return
     */
    public double distTo(Point destination) {
      if (!hasPathTo(destination))
        throw new IllegalArgumentException("no path");
      return distTo[destination.getX()][destination.getY()];
    }
    public boolean hasPathTo(Point destination) {
      return distTo[destination.getX()][destination.getY()] != Double.POSITIVE_INFINITY;
    }
    public Collection<Point> pathTo(Point destination) {
      if (!hasPathTo(destination))
        throw new IllegalArgumentException("no path");
      Deque<Point> path = new ArrayDeque<>();
      for (Point p = destination; p != null; p = closestParentOf[p.getX()][p.getY()]) {
        path.addFirst(p);
      }
      return path;
    }
    private void relax(Point parent, Point child) {
      if (distTo[child.getX()][child.getY()] > distTo[parent.getX()][parent.getY()] + matrix[child.getX()][child.getY()]) {
        distTo[child.getX()][child.getY()] = distTo[parent.getX()][parent.getY()] + matrix[child.getX()][child.getY()];
        closestParentOf[child.getX()][child.getY()] = parent;
      }
    }
  }

  private interface NeighborFinder {
    /**
     * return neighbors of point(x,y)
     * @param p
     * @return
     */
    Iterable<Point> adj(Point p);
  }
  private class Point {
    private final int x;
    private final int y;
    public int getX() {
      return x;
    }

    public int getY() {
      return y;
    }

    public Point(int x, int y) {
      this.x = x;
      this.y = y;
    }
    @Override
    public boolean equals(Object o) {
      if (o == null)
        return false;
      if (o.getClass() != Point.class)
        return false;
      Point that = (Point) o;
      return this.x == that.x && this.y == that.y;
    }
    @Override
    public int hashCode() {
      int result = 17;
      result = result * 31 + x;
      result = result * 31 + y;
      return result;
    }
  }
  private class VerticalNeighborFinder implements NeighborFinder{
    private int width;
    private int height;
    public VerticalNeighborFinder(int height, int width) {
      this.width = width;
      this.height = height;
    }
    @Override
    public Iterable<Point> adj(Point p) {
      int x = p.getX();
      int y = p.getY();
      List<Point> neighbors = new ArrayList<>();
      if (y == height - 1)
        return neighbors;
      if (width == 1) {
        neighbors.add(new Point(x, y + 1));
      } else if (x == 0) {
        neighbors.add(new Point(x, y + 1));
        neighbors.add(new Point(x + 1, y + 1));
      } else if (x == width - 1) {
        neighbors.add(new Point(x - 1, y + 1));
        neighbors.add(new Point(x, y + 1));
      } else {
        neighbors.add(new Point(x - 1, y + 1));
        neighbors.add(new Point(x, y + 1));
        neighbors.add(new Point(x + 1, y + 1));
      }
      return neighbors;
    }
  }
  private class HorizontalNeighborFinder implements NeighborFinder{
    private int width;
    private int height;
    public HorizontalNeighborFinder(int height, int width) {
      this.width = width;
      this.height = height;
    }
    @Override
    public Iterable<Point> adj(Point p) {
      int x = p.getX();
      int y = p.getY();
      List<Point> neighbors = new ArrayList<>();
      if (x == width - 1)
        return neighbors;
      if (height == 1) {
        neighbors.add(new Point(x + 1, y));
      } else if (y == 0) {
        neighbors.add(new Point(x + 1, y));
        neighbors.add(new Point(x + 1, y + 1));
      } else if (y == height - 1) {
        neighbors.add(new Point(x + 1, y));
        neighbors.add(new Point(x + 1, y - 1));
      } else {
        neighbors.add(new Point(x + 1, y - 1));
        neighbors.add(new Point(x + 1, y));
        neighbors.add(new Point(x + 1, y + 1));
      }
      return neighbors;
    }
  }
}
