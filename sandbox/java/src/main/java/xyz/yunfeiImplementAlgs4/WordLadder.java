package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 2/4/17.
 */

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.SET;
import edu.princeton.cs.algs4.StdIn;

import java.util.*;

/******************************************************************************
 *  Compilation:  javac WordLadder.java
 *  Execution:    java WordLadder word1 word2 < wordlist.txt
 *  Dependencies: Graph.java IndexSET.java In.java BreadthFirstPaths.java
 *
 *  Data files:   http://algs4.cs.princeton.edu/41graph/words5.txt
 *                http://algs4.cs.princeton.edu/41graph/words6.txt
 *                http://algs4.cs.princeton.edu/41graph/words5-knuth.txt
 *
 *  Creates a minimum length word ladder connecting two words.
 *
 *  java WordLadder words5.txt
 *  flirt break
 *  length = 11
 *  flirt
 *  flint
 *  fling
 *  cling
 *  clink
 *  click
 *  clock
 *  cloak
 *  croak
 *  creak
 *  break
 *
 *  allow brown
 *  NOT CONNECTED
 *
 *  white house
 *  length = 18
 *  white
 *  while
 *  whale
 *  shale
 *  shake
 *  slake
 *  slate
 *  plate
 *  place
 *  peace
 *  peach
 *  poach
 *  coach
 *  couch
 *  cough
 *  rough
 *  rouge
 *  rouse
 *  house
 *
 *  % java WordLadder words5-knuth.txt
 *  white house
 *  length = 9
 *  white
 *  whits
 *  shits
 *  shots
 *  soots
 *  roots
 *  routs
 *  route
 *  rouse
 *  house
 *
 ******************************************************************************/
public class WordLadder {
  private Graph g;
  private Map<String, Integer> word2Index;
  private List<String> index2Word;
  private static Iterable<String> emptyPath = Arrays.asList(new String[0]);


  /**
   * add a word to the ladder
   * @param w
   */
  private void addWord(String w) {
    if (word2Index.containsKey(w)) {
      return;
    }
    word2Index.put(w, word2Index.size());
    index2Word.add(w);
    char[] letters = w.toCharArray();
    for (int i = 0; i < w.length(); i++) {
      char oldChar = letters[i];
      assert(oldChar <= 'z' && oldChar >= 'a');
      for (char c = 'a'; c <= 'z'; c++) {
        if (oldChar == c) {
          continue;
        }
        letters[i] = c;
        if (word2Index.containsKey(new String(letters)))
          g.addEdge(word2Index.get(w), word2Index.get(new String(letters)));
      }
      letters[i] = oldChar;
    }
  }
  /**
   *
   */
  public WordLadder(In in) {
    Set<String> words = new HashSet<String>();
    while (!in.isEmpty()) {
      String word = in.readString();
      words.add(word);
    }
    System.err.println("Finished reading word list");
    g = new Graph(words.size());
    word2Index = new HashMap<>();
    index2Word = new ArrayList<>();
    for (String w : words)
      addWord(w);
  }

  public Iterable<String> getPath(String from, String to) {
    if (!word2Index.containsKey(from) || !word2Index.containsKey(to)) {
      return Arrays.asList(new String[0]);
    }
    BreadthFirstPaths dfs = new BreadthFirstPaths(g, word2Index.get(from));
    List<String> path = new ArrayList<>();
    for (int i : dfs.pathTo(word2Index.get(to))) {
      path.add(index2Word.get(i));
    }
    return path;
  }

  public static void main(String[] args) {
    In in = new In(args[0]);
    WordLadder ladder = new WordLadder(in);
    String from = args[1];
    String to = args[2];
    for (String w : ladder.getPath(from, to)) {
      System.out.print(w + " ");
    }
    System.out.print("\n");
  }
}
