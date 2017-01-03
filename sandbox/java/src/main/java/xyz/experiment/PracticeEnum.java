package experiment;

import java.util.EnumSet;
import java.util.Set;

/**
 * Created by guoy28 on 11/2/16.
 */
public class PracticeEnum {
  private enum Color {
    RED(255, 0, 0), GREEN(0, 255, 0), BLUE(0, 0, 255);
    private int r, g, b;
    private Color(int r, int g, int b) {
      this.r = r;
      this.g = g;
      this.b = b;
    }
    public int getR() {
      return r;
    }
    public int getG() {
      return g;
    }
    public int getB() {
      return b;
    }
  }

  public static void main(String[] args) {
    EnumSet<Color> pink = EnumSet.of(Color.RED, Color.BLUE);
    EnumSet<Color> white = EnumSet.of(Color.RED, Color.BLUE, Color.GREEN);
    EnumSet<Color> yellow = EnumSet.of(Color.BLUE, Color.GREEN);
    draw(pink);
    draw(white);
    draw(yellow);
  }
  public static void draw(Set<Color> compositeColor) {
    for (Color c : compositeColor) {
      System.out.print(c + " ");
    }
    System.out.println();
  }
}
