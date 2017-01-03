package experiment;

import java.util.StringTokenizer;

/**
 * Created by guoy28 on 12/18/16.
 */
public class Tokenizer {
  public static void main(String[] args) {
    String s = "1\t2\t\t4";
    StringTokenizer tok = new StringTokenizer(s);
    while(tok.hasMoreTokens()) {
      System.out.println(tok.nextToken());
    }
  }
}
