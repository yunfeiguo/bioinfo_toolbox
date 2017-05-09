package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/19/17.
 */
/******************************************************************************
 *  Compilation:  javac Genome.java
 *  Execution:    java Genome - < input.txt   (compress)
 *  Execution:    java Genome + < input.txt   (expand)
 *  Dependencies: BinaryIn.java BinaryOut.java
 *  Data files:   http://algs4.cs.princeton.edu/55compression/genomeTiny.txt
 *
 *  Compress or expand a genomic sequence using a 2-bit code.
 *
 *  % more genomeTiny.txt
 *  ATAGATGCATAGCGCATAGCTAGATGTGCTAGC
 *
 *  % java Genome - < genomeTiny.txt | java Genome +
 *  ATAGATGCATAGCGCATAGCTAGATGTGCTAGC
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.*;

/**
 *  The {@code Genome} class provides static methods for compressing
 *  and expanding a genomic sequence using a 2-bit code.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/55compress">Section 5.5</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Genome {
  private Genome() {}
  public static void compress(BinaryIn in, BinaryOut out) {
    while(!in.isEmpty()) {
      char c = in.readChar();
      try {
        int i = edu.princeton.cs.algs4.Alphabet.DNA.toIndex(c);
        out.write(i, 2);
      } catch (IllegalArgumentException e) {
        continue;
      }
    }
  }
  public static void expand(BinaryIn in, BinaryOut out) {
    while(!in.isEmpty()) {
      int i = in.readInt(2);
      char c = edu.princeton.cs.algs4.Alphabet.DNA.toChar(i);
      out.write(c);
    }
  }
  /**
   * Sample client that calls {@code compress()} if the command-line
   * argument is "-" an {@code expand()} if it is "+".
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    BinaryIn in = new BinaryIn(args[1]);
    BinaryOut out = new BinaryOut(args[2]);
    if      (args[0].equals("-")) {
      compress(in, out);
    } else if (args[0].equals("+")) {
      expand(in, out);
    } else throw new IllegalArgumentException("Illegal command line argument");
    out.close();
  }
}
