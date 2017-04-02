package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/2/17.
 */
/******************************************************************************
 *  Compilation:  javac RabinKarp.java
 *  Execution:    java RabinKarp pat txt
 *  Dependencies: StdOut.java
 *
 *  Reads in two strings, the pattern and the input text, and
 *  searches for the pattern in the input text using the
 *  Las Vegas version of the Rabin-Karp algorithm.
 *
 *  % java RabinKarp abracadabra abacadabrabracabracadabrabrabracad
 *  pattern: abracadabra
 *  text:    abacadabrabracabracadabrabrabracad
 *  match:                 abracadabra
 *
 *  % java RabinKarp rab abacadabrabracabracadabrabrabracad
 *  pattern: rab
 *  text:    abacadabrabracabracadabrabrabracad
 *  match:           rab
 *
 *  % java RabinKarp bcara abacadabrabracabracadabrabrabracad
 *  pattern: bcara
 *  text:         abacadabrabracabracadabrabrabracad
 *
 *  %  java RabinKarp rabrabracad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern:                        rabrabracad
 *
 *  % java RabinKarp abacad abacadabrabracabracadabrabrabracad
 *  text:    abacadabrabracabracadabrabrabracad
 *  pattern: abacad
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.StdOut;

/**
 *  The {@code RabinKarp} class finds the first occurrence of a pattern string
 *  in a text string.
 *  <p>
 *  This implementation uses the Rabin-Karp algorithm.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/53substring">Section 5.3</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 */

public class RabinKarp {
  private static final int prime = 997;
  private static final int radix = 3;
  private int patternHash;
  private String pattern;
  public RabinKarp(String pattern) {
    patternHash = hash(pattern);
    this.pattern = pattern;
  }
  private int hash(String s) {
    int hash = 0;
    /*
    hash = (s_0*radix^(L-1) + s_1*radix^(L-2) + ... + s_(L-1)) % prime
         = (((s_0*radix + s_1)*radix + s_2)*radix ...) % prime
    (ab)%prime = ((a%prime) * (b%prime))%prime
    (a + b) % prime = (a%prime + b%prime)%prime
     */
    for (int i = 0; i < s.length(); i++) {
      hash = ((hash * (radix % prime))%prime + s.charAt(i)) % prime;
    }
    return hash;
  }

  /**
   *
   * @param base
   * @param power
   * @param mod
   * @return (base^power)%mod
   */
  private int powerMod(int base, int power, int mod) {
    int result = 1;
    for (int i = 0; i < power; i++) {
      result = (result * (base % mod)) % mod;
    }
    return result;
  }

  /**
   * Rabin-Karp algorithm to find first occurence of pattern in txt
   * @param txt
   * @return
   */
  public int search(String txt) {
    if (txt.length() < pattern.length()) throw new IllegalArgumentException("query string too short");
    int currentHash = hash(txt.substring(0, pattern.length()));
    if (currentHash == patternHash && pattern.compareTo(txt.substring(0, pattern.length())) == 0)
      return 0;
    int multiplier = powerMod(radix, pattern.length() - 1, prime);
    /*
    assume we know
    x_i = (txt_i*radix^(L-1) + txt_(i+1)*radix^(L-2) + ... + txt_(i+L-1)) % prime
    we want to calculate x_(i+1) in constant time

    x_(i+1) = (txt_(i+1)*radix^(L-1) + txt_(i+2)*radix^(L-2) + ... + txt_(i+L)) % prime
            = (((x_i - txt_i*radix^(L-1)%prime)*(radix%prime))%prime + txt_(i+L)) % prime
     */
    for (int i = 1; i < txt.length() - pattern.length() + 1; i++) {
      currentHash = (
              ((currentHash - (multiplier * (txt.charAt(i-1)%prime))%prime)*
                      (radix%prime))%prime +
              txt.charAt(i + pattern.length() - 1)
              ) % prime;
      /*
      java % can operate on negative numbers, but
      we assume positive numbers all the time
       */
      if (currentHash < 0)
        currentHash += prime;
      if (currentHash == patternHash && pattern.compareTo(txt.substring(i, i + pattern.length())) == 0)
        return i;
    }
    return txt.length();
  }

  /**
   * Takes a pattern string and an input string as command-line arguments;
   * searches for the pattern string in the text string; and prints
   * the first occurrence of the pattern string in the text string.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    String pat = args[0];
    String txt = args[1];

    RabinKarp searcher = new RabinKarp(pat);
    int offset = searcher.search(txt);

    // print results
    StdOut.println("text:    " + txt);

    // from brute force search method 1
    StdOut.print("pattern: ");
    for (int i = 0; i < offset; i++)
      StdOut.print(" ");
    StdOut.println(pat);
  }
}
