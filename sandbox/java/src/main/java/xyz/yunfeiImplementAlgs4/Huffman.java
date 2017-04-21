package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/19/17.
 */
/******************************************************************************
 *  Compilation:  javac Huffman.java
 *  Execution:    java Huffman - < input.txt   (compress)
 *  Execution:    java Huffman + < input.txt   (expand)
 *  Dependencies: BinaryIn.java BinaryOut.java
 *  Data files:   http://algs4.cs.princeton.edu/55compression/abra.txt
 *                http://algs4.cs.princeton.edu/55compression/tinytinyTale.txt
 *                http://algs4.cs.princeton.edu/55compression/medTale.txt
 *                http://algs4.cs.princeton.edu/55compression/tale.txt
 *
 *  Compress or expand a binary input stream using the Huffman algorithm.
 *
 *  % java Huffman - < abra.txt | java BinaryDump 60
 *  010100000100101000100010010000110100001101010100101010000100
 *  000000000000000000000000000110001111100101101000111110010100
 *  120 bits
 *
 *  % java Huffman - < abra.txt | java Huffman +
 *  ABRACADABRA!
 *
 ******************************************************************************/

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;

import java.util.*;

/**
 *  The {@code Huffman} class provides static methods for compressing
 *  and expanding a binary input using Huffman codes over the 8-bit extended
 *  ASCII alphabet.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/55compress">Section 5.5</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Huffman {
  private static final int R = 256; //alphabet size
  private static final int W = 8; //# of bits for char
  private static final char NULL = '\0'; //terminal char
  private Huffman() { }
  private static class Node implements Comparable<Node>{
    public Node(Node left, Node right, char c, int weight) {
      this.left = left;
      this.right = right;
      this.c = c;
      this.weight = weight;
    }
    public Node left;
    public Node right;
    public char c;
    public int weight;

    @Override
    public int compareTo(Node o) {
      return this.weight - o.weight;
    }
  }
  public static void compress(BinaryIn in, BinaryOut out) {
    //not space efficient, good for now
    List<Character> input = new ArrayList<>();
    Node root = constructHuffmanTrie(in, input);
    writeHuffmanTrie(root, out);
    //save huffman codes in a symbol table
    Map<Character, Boolean[]> st = new HashMap<>();
    Deque<Boolean> codeOnStack = new ArrayDeque<>();
    //in case there is only ONE unique char in input
    if (root.left == null && root.right == null) {
      st.put(root.c, new Boolean[]{false});
    } else {
      getHuffmanCode(root, st, codeOnStack);
    }
    //write # of chars, why?
    //because we use variable length codes, we don't when we
    //end, there could be padding at the end
    out.write(input.size());
    //write huffman code
    for (char c : input) {
      for (boolean b : st.get(c)) {
        out.write(b);
      }
    }
  }

  /**
   * construct Huffman trie and store input chars
   * @param in
   * @param input
   */
  private static Node constructHuffmanTrie(BinaryIn in, List<Character> input) {
    PriorityQueue<Node> q = new PriorityQueue<>();
    //not space efficient, good for now
    int[] frequency = new int[R];
    while(!in.isEmpty()) {
      char c = in.readChar(W);
      frequency[c]++;
      input.add(c);
    }
    for (int i = 0; i < frequency.length; i++) {
      if (frequency[i] > 0)
        q.offer(new Node(null, null, (char) i, frequency[i]));
    }
    //merging nodes from lowest weight to highest
    while(q.size() > 1) {
      Node left = q.poll();
      Node right = q.poll();
      Node parent = new Node(left, right, NULL, left.weight + right.weight);
      q.offer(parent);
    }
    return q.isEmpty()? null : q.poll();
  }

  /**
   * store huffman code in a symbol table
   * as array of booleans
   * @param n
   * @param st
   * @param stack
   */
  private static void getHuffmanCode(Node n, Map<Character, Boolean[]> st, Deque<Boolean> stack) {
    if (n.left == null && n.right == null) {
      st.put(n.c, stack.toArray(new Boolean[0]));
      return;
    }
    if (n.left != null) {
      stack.addLast(false);
      getHuffmanCode(n.left, st, stack);
      stack.removeLast();
    }
    if (n.right != null) {
      stack.addLast(true);
      getHuffmanCode(n.right, st, stack);
      stack.removeLast();
    }
  }
  /**
   * implicitly assume all leaf nodes have a char associated with it
   */
  private static void writeHuffmanTrie(Node n, BinaryOut out) {
    if (n.left == null && n.right == null) {
      out.write(true); //mark a leaf node
      out.write(n.c, W);
      return;
    }
    out.write(false);
    if (n.left != null) {
      writeHuffmanTrie(n.left, out);
    }
    if (n.right != null) {
      writeHuffmanTrie(n.right, out);
    }
    if ((n.left == null) ^ (n.right == null))
      throw new RuntimeException("child nodes must both be null or not");
  }
  public static void expand(BinaryIn in, BinaryOut out) {
    Node root = readHuffmanTrie(in);
    int l = in.readInt();
    int count = 0;
    Node current = root;
    while(!in.isEmpty() && count < l) {
      boolean b = in.readBoolean();
      if (current.left == null && current.right == null) {
        out.write(current.c, W);
        count++;
        current = root;
      }
      current = b ? current.right :current.left;
      if (current == null)
        //in case there is only one char in decompressed text
        current = root;
    }
  }
  private static Node readHuffmanTrie(BinaryIn in) {
    boolean b = in.readBoolean();
    if (b) {
      return new Node(null, null, in.readChar(W), 0);
    }
    return new Node(readHuffmanTrie(in), readHuffmanTrie(in), NULL, 0);
  }
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
