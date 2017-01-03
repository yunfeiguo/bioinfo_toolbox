package algorithm_practice.leetcode;

/**
 * Created by guoy28 on 12/8/16.
 */
public class SumOfTwo {
  public static void main(String[] args) {
    System.out.println(getSum(Integer.MIN_VALUE,200) == Integer.MIN_VALUE + 200);
    System.out.println(getSum(Integer.MAX_VALUE,Integer.MIN_VALUE) == Integer.MAX_VALUE + Integer.MIN_VALUE);
    System.out.println(getSum(-99, 200) == -99 + 200);
    System.out.println(getSum(-99, -9) == -99 + -9);
  }
  public static int getSum(int a, int b) {
    //assume no overflow
    int advance = 0;
    int result = 0;
    int aFlip =  a & ((1<<31)-1);
    int bFlip =  b & ((1<<31)-1);
    System.out.println("a: " + Integer.toBinaryString(a));
    System.out.println("a leftmost 0: " + aFlip);
    System.out.println("a leftmost 0: " + Integer.toBinaryString(aFlip));
    System.out.println("b: " + Integer.toBinaryString(b));
    System.out.println("b leftmost 0: " + bFlip);
    System.out.println("b leftmost 0: " + Integer.toBinaryString(bFlip));
    //from right to left
    for (int i = 0; i < 32; i++) {
      result |= (getBit(a, i) ^ getBit(b, i) ^ advance) << i;
      advance = (getBit(a, i) & getBit(b, i)) | (getBit(a, i) & advance) | (getBit(b, i) & advance);
    }
    System.out.println("result: " + Integer.toBinaryString(result));
    System.out.println("result: " + result);
    return result;
  }
  //get ith bit in n (from right)
  public static int getBit(int n, int i) {
    return ((1 << i) & n) == 0? 0 : 1;
  }
}
