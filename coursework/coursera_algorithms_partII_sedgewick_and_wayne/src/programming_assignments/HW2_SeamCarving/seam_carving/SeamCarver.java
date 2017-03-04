package programming_assignments.HW2_SeamCarving.seam_carving;

import com.sun.tools.corba.se.idl.constExpr.Not;
import edu.princeton.cs.algs4.Picture;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

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
    pic = picture;
    energy = new double[width()][height()];
    for (int x = 0; x < width(); x++) {
      for (int y = 0; y < height(); y++) {
        energy[x][y] = energy(x, y);
      }
    }
  }
  // current picture
  public Picture picture() {
    return pic;
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
      return(BORDER_ENERGY);
    Color left = pic.get(x - 1, y);
    Color right = pic.get(x - 1, y);
    Color top = pic.get(x - 1, y);
    Color bottom = pic.get(x - 1, y);
    double rX = left.getRed() - right.getRed();
    double gX = left.getGreen() - right.getGreen();
    double bX = left.getBlue() - right.getBlue();
    double rY = top.getRed() - bottom.getRed();
    double gY = top.getGreen() - bottom.getGreen();
    double bY = top.getBlue() - bottom.getBlue();
    double deltaX = rX * rX + gX * gX + bX * bX;
    double deltaY = rY * rY + gY * gY + bY * bY;
    return(Math.sqrt(deltaX * deltaX + deltaY + deltaY));
  }

  /**
   *
   * sequence of indices for horizontal seam
   * seam has minimum energy among all possible
   * seams
    */
  public int[] findHorizontalSeam() {
    NeighborFinder neighborFinder = new HorizontalNeighborFinder(height(), width());
    double minEnergy = Double.POSITIVE_INFINITY;
    int minEnergyColumn;
    List<Point> firstRow = new ArrayList<>(width());
    for (int x = 0; x < width(); x++) {
      firstRow.add(new Point(x, 0));
    }
    MatrixDigraphSP sp = new MatrixDigraphSP(energy, neighborFinder, firstRow);
    for (int x = 0; x < width(); x++) {
      if (sp.distTo(new Point(x, height() - 1)) < minEnergy) {
        minEnergy = sp.distTo(new Point(x, height() - 1));
        minEnergyColumn = x;
      }
    }
    return(sp.pathTo(new Point(minEnergyColumn, height() - 1)));
  }

  /** sequence of indices for vertical seam
   * seam has minimum energy among all possible
   * seams
   */
  public int[] findVerticalSeam() {
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
      if (i < 0 || i >= width())
        throw new IllegalArgumentException("seam index out of bound");
    }
    for (int i = 1; i < seam.length - 1; i++) {
      if (Math.abs(seam[i] - seam[i - 1]) > 1 || Math.abs(seam[i + 1] - seam[i]) > 1)
        throw new IllegalArgumentException("neighbors not differ <= 1");
    }
    if (Math.abs(seam[1] - seam[0]) > 1 || Math.abs(seam[width() - 1] - seam[width() - 2]) > 1)
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
    pic = newPic;
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
      if (i < 0 || i >= height())
        throw new IllegalArgumentException("seam index out of bound");
    }
    for (int i = 1; i < seam.length - 1; i++) {
      if (Math.abs(seam[i] - seam[i - 1]) > 1 || Math.abs(seam[i + 1] - seam[i]) > 1)
        throw new IllegalArgumentException("neighbors not differ <= 1");
    }
    if (Math.abs(seam[1] - seam[0]) > 1 || Math.abs(seam[height() - 1] - seam[height() - 2]) > 1)
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
    pic = newPic;
  }

  public static void main(String[] args) {
  }
  /**
   * matrix-based directed graph
   * edges defined dynamically
   * find shortest path (SP) from sources
   * to all possible destinations
   */
  private class MatrixDigraphSP {
    private double[][] matrix;
    private List<List<Double>> distTo;
    private List<List<Point>> parentOf;

    /**
     * constructor will find shortest paths from sources
     * to all possible destinations
     * @param matrix
     * @param neighborFinder
     * @param sources
     */
    public MatrixDigraphSP(double[][] matrix, NeighborFinder neighborFinder, Iterable<Point> sources) {

    }

    /**
     * return
     * @return
     */
    public double distTo(Point destination) {
      if (!hasPathTo(destination))
        throw new IllegalArgumentException("no path");
      return (double) -1;
    }
    public boolean hasPathTo(Point destination) {
      throw new NotImplementedException();
      return distTo.get(destination.getX()).get(destination.getY()) != Double.POSITIVE_INFINITY;
    }
    public int[] pathTo(Point destination) {
      if (!hasPathTo(destination))
        throw new IllegalArgumentException("no path");
      return new int[0];
    }
  }

  private interface NeighborFinder {
    /**
     * return neighbors of point(x,y)
     * @param x
     * @param y
     * @return
     */
    public Iterable<Point> adj(int x, int y);
  }
  private class Point {
    public int getX() {
      return x;
    }

    public int getY() {
      return y;
    }
    private final int x;
    private final int y;
    public Point(int x, int y) {
      this.x = x;
      this.y = y;
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
    public Iterable<Point> adj(int x, int y) {
      throw new NotImplementedException();
      List<Point> neighbors = new ArrayList<>();
      if (x == width - 1)
        return neighbors;
      if (y == 0) {
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
  private class HorizontalNeighborFinder implements NeighborFinder{
    private int width;
    private int height;
    public HorizontalNeighborFinder(int height, int width) {
      this.width = width;
      this.height = height;
    }
    @Override
    public Iterable<Point> adj(int x, int y) {
      List<Point> neighbors = new ArrayList<>();
      if (x == width - 1)
        return neighbors;
      if (y == 0) {
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
