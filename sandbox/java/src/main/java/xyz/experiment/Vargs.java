package main.java.xyz.experiment;

/**
 * Created by guoy28 on 1/13/17.
 */
public class Vargs {
  public static int sum(int... x) {
    int sum = 0;
    for (int i : x) {
      sum += i;
    }
    return sum;
  }

  public static void main(String[] args) {
    System.out.println(sum(1,2,3,4));
  }
}
