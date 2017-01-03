package algorithm_practice.leetcode;

/**
 * Created by guoy28 on 12/2/16.
 */
public class NimGame {
    public boolean canWinNim(int n) {
      //assume n>=1
      if (n <= 3) {
        return true;
      }
      boolean[] canWin = new boolean[3];
      canWin[0] = true;
      canWin[1] = true;
      canWin[2] = true;
      for (int i = 3; i < n; i++) {
        canWin[i % 3] = !canWin[0] || !canWin[1] || !canWin[2];
      }
      return canWin[(n-1) % 3];
    }

  public static void main(String[] args) {
    NimGame test = new NimGame();
    System.out.println(!test.canWinNim(1348820612));
    System.out.println((Integer.MAX_VALUE) * ((double) 1 /Long.MAX_VALUE));
  }
}
