package algorithm_practice.leetcode;

/**
 * Created by guoy28 on 10/30/16.
 */
public class WildcardMatching {
    public boolean isMatch(String s, String p) {
      if (s == null || p == null) {
        return false;
      }
      char[] x = s.toCharArray();
      char[] y = p.toCharArray();
      return isMatch(x, y, 0, 0);
    }
    private boolean isMatch(char[] s, char[] p, int sStart, int pStart) {
      if (sStart == s.length && pStart == p.length) {
        return true;
      }
      if ((sStart < s.length && pStart < p.length) && (s[sStart] == p[pStart] || p[pStart] == '?')) {
        if (isMatch(s, p, sStart + 1, pStart + 1)) {
          return true;
        }
      } else if (pStart < p.length && p[pStart] == '*') {
        if (isMatch(s, p, sStart, pStart + 1)) {
          return true;
        }
        if (sStart < s.length) {
          if (isMatch(s, p, sStart + 1, pStart + 1) || isMatch(s, p, sStart + 1, pStart)) {
            return true;
          }
        }
      }
      return false;
    }

  public static void main(String[] args) {
    WildcardMatching test = new WildcardMatching();
    System.out.println(test.isMatch(
            "babaaababaabababbbbbbaabaabbabababbaababbaaabbbaaab",
            "***bba**a*bbba**aab**b") == false);
  }
}
