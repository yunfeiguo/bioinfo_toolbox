package programming_assignments.HW2_SeamCarving.hw2_SeamCarving;

import edu.princeton.cs.algs4.Picture;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.awt.*;

/**
 * Created by guoy28 on 3/4/17.
 */
public class SeamCarver {
  private static int BORDER_ENERGY = 1000;
  private int height;
  private int width;
  private Picture pic;
  // create a seam carver object based on the given picture
  public SeamCarver(Picture picture) {
    pic = picture;
  }
  // current picture
  public Picture picture() {
    return pic;
  }
  // width of current picture
  public int width() {
    return width;
  }
  // height of current picture
  public int height() {
    return height;
  }
  // energy of pixel at column x and row y
  public double energy(int x, int y) {
    if (x < 0 || x >= width || y < 0 || y >= height)
      throw new IndexOutOfBoundsException("coordinates out of bound");
    if (x == 0 || x == width - 1 || y == 0 || y == height - 1)
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
  // sequence of indices for horizontal seam
  public int[] findHorizontalSeam() {
  }
  // sequence of indices for vertical seam
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
    if (seam.length != width)
      throw new IllegalArgumentException("unequal to width");
    if (height <= 1)
      throw new IllegalArgumentException("image too small");
    for (int i : seam) {
      if (i < 0 || i >= width)
        throw new IllegalArgumentException("seam index out of bound");
    }
    for (int i = 1; i < seam.length; i++) {
      throw new NotImplementedException();
    }
  }

  /** remove vertical seam from current picture
   *
   * @param seam
   */
  public void removeVerticalSeam(int[] seam) {
    if (seam == null)
      throw new NullPointerException("nothing to remove");
    if (seam.length != height)
      throw new IllegalArgumentException("unequal to height");
    if (width <= 1)
      throw new IllegalArgumentException("image too small");
    for (int i : seam) {
      if (i < 0 || i >= height)
        throw new IllegalArgumentException("seam index out of bound");
    }
    for (int i = 1; i < seam.length; i++) {
      throw new NotImplementedException();
    }
  }
}
