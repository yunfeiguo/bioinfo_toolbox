package programming_assignments.HW4_boggle.boggle;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.TST;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by guoy28 on 4/5/17.
 */
public class BoggleSolver {
  private Node root;
  /**
  // Initializes the data structure using the given array of strings as the dictionary.
  // (You can assume each word in the dictionary contains only the uppercase letters A through Z.)
   */
  public BoggleSolver(String[] dictionary) {
    /*
    construct a ternary trie to store the strings in dict
     */
    root = null;
    for (String s : dictionary) {
      root = insert(s, root, 0);
    }
  }

  /**
   * consider n as root, insert
   * string s into ternary trie
   * based on s[i] character
   *
   * @param s
   * @param n
   * @param i
   * @return
   */
  private Node insert(String s, Node n, int i) {
    if (n == null) {
      n = new Node();
      n.key = s.charAt(i);
    }
    if (i == s.length() - 1 && n.key == s.charAt(i)) {
      n.value = s;
      return n;
    }
    if (s.charAt(i) < n.key) {
      n.left = insert(s, n.left, i);
    } else if (s.charAt(i) > n.key) {
      n.right = insert(s, n.right, i);
    } else {
      n.middle = insert(s, n.middle, i + 1);
    }
    return n;
  }

  /**
   * start from n, find s[i:]
   *
   * @param s
   * @param n
   * @param i
   * @return
   */
  private Node find(String s, Node n, int i) {
    if (s == null || n == null || i >= s.length())
      return null;
    if (s.charAt(i) > n.key)
      return find(s, n.right, i);
    else if (s.charAt(i) < n.key)
      return find(s, n.left, i);
    else if (i == s.length() - 1)
      return s.compareTo(n.value) == 0? n : null;
    else
      return find(s, n.middle, i+1);
  }
  private class Node {
    public char key;
    public String value;
    public Node left;
    public Node middle;
    public Node right;
  }

  /**
  // Returns the set of all valid words in the given Boggle board, as an Iterable.
   */
  public Iterable<String> getAllValidWords(BoggleBoard board) {
    int nrow = board.rows();
    int ncol = board.cols();
    boolean[][] visited = new boolean[nrow][ncol];
    Set<String> validWords = new HashSet<>();

    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        //do we need to check whether duplicates exist?
        getAllValidWordsFrom(i, j, root, visited, validWords, board);
      }
    }
    return validWords;
  }
  private void getAllValidWordsFrom(int i, int j, Node n, boolean[][] visited, Set<String> validWords, BoggleBoard b) {
    if (i < 0 || i >= b.rows() || j < 0 || j >= b.cols() || visited[i][j] || n == null)
      return;
    char currentChar = b.getLetter(i, j);
    if (currentChar > n.key) {
      getAllValidWordsFrom(i, j, n.right, visited, validWords, b);
    } else if (currentChar < n.key) {
      getAllValidWordsFrom(i, j, n.left, visited, validWords, b);
    } else {
      visited[i][j] = true;
      if (n.value != null)
        validWords.add(n.value);
      n = n.middle;
      //QU always together
      if (currentChar == 'Q') {
        if (n.key == 'U')
          n = n.middle;
        else {
          visited[i][j] = false;
          return;
        }
      }
      /*
      ↖ ↑ ↗
       ← · →
      ↙ ↓ ↘
      8 directions allowed
       */
      getAllValidWordsFrom(i - 1, j, n, visited, validWords, b);
      getAllValidWordsFrom(i + 1, j, n, visited, validWords, b);
      getAllValidWordsFrom(i, j - 1, n, visited, validWords, b);
      getAllValidWordsFrom(i, j + 1, n, visited, validWords, b);

      getAllValidWordsFrom(i - 1, j - 1, n, visited, validWords, b);
      getAllValidWordsFrom(i - 1, j + 1, n, visited, validWords, b);
      getAllValidWordsFrom(i + 1, j - 1, n, visited, validWords, b);
      getAllValidWordsFrom(i + 1, j + 1, n, visited, validWords, b);
      visited[i][j] = false;
    }
    return;
  }

  /**
  // Returns the score of the given word if it is in the dictionary, zero otherwise.
  // (You can assume the word contains only the uppercase letters A through Z.)
  */
  public int scoreOf(String word) {
    if (find(word, root, 0) == null)
      return 0;
    int l = word.length();
    if (l <=2 )
      return 0;
    else if (l >= 3 && l <= 4)
      return 1;
    else if (l == 5)
      return 2;
    else if (l == 6)
      return 3;
    else if (l == 7)
      return 5;
    else if (l >= 8)
      return 11;
    return 0;
  }
  public static void main(String[] args) {
    In in = new In(args[0]);
    String[] dictionary = in.readAllStrings();
    BoggleSolver solver = new BoggleSolver(dictionary);
    BoggleBoard board = new BoggleBoard(args[1]);
    int score = 0;
    for (String word : solver.getAllValidWords(board)) {
      StdOut.println(word);
      score += solver.scoreOf(word);
    }
    StdOut.println("Score = " + score);
    //benchmark
    long start = System.nanoTime();
    int n = 10000;
    for (int i = 0; i < n; i++) {
      solver.getAllValidWords(board);
    }
    long end = System.nanoTime();
    StdOut.println("Time for " + n + " searches: " + (end - start)/1e9 + " seconds");
  }
}
