package algorithm_practice.leetcode;

/**
 * Created by guoy28 on 10/26/16.
 */

import javax.management.RuntimeErrorException;
import java.util.Random;

/**
 * suppose we have a method called rand5 that can be used for
 * generating random numbers from 0~4, how to get a method
 * for generating random numbers from 0~6?
 *
 */
public class Rand7FromRand5 {
  private Random randomNumberGenerator;

  public Rand7FromRand5(int seed) {
    this.randomNumberGenerator = new Random(seed);
  }

  public int rand7() {
    for(;;) {
      int next = rand5() * 5 + rand5();
      if (next >= 3*7) {
        continue;
      }
      return next % 7;
    }
  }

  public int rand5() {
    return randomNumberGenerator.nextInt(5);
  }

  public static void main(String[] args) {
    Rand7FromRand5 test = new Rand7FromRand5(11);
    int total = 70000;
    int[] count = new int[7];
    //simulate some random numbers
    for (int i = 0; i < total;i++) {
      count[test.rand7()]++;
    }
    //check if they are truly random
    double tolerance = 0.01;
    for (int n : count) {
      if (Math.abs((double) n / total - 1.0/7) > tolerance) {
        throw new RuntimeException("frequency deviates from expectation!");
      }
    }
  }
}
