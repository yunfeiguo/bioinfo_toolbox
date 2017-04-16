package yunfeiImplementAlgs4;

import java.util.*;

import edu.princeton.cs.algs4.*;
import edu.princeton.cs.algs4.DepthFirstSearch;
/******************************************************************************
 *  Compilation:  javac NFA.java
 *  Execution:    java NFA regexp text
 *  Dependencies: Stack.java Bag.java Digraph.java DirectedDFS.java
 *
 *  % java NFA "(A*B|AC)D" AAAABD
 *  true
 *
 *  % java NFA "(A*B|AC)D" AAAAC
 *  false
 *
 *  % java NFA "(a|(bc)*d)*" abcbcd
 *  true
 *
 *  % java NFA "(a|(bc)*d)*" abcbcbcdaaaabcbcdaaaddd
 *  true
 *
 *  Remarks
 *  -----------
 *  The following features are not supported:
 *    - The + operator
 *    - Multiway or
 *    - Metacharacters in the text
 *    - Character classes.
 *
 ******************************************************************************/

/**
 *  The {@code NFA} class provides a data type for creating a
 *  <em>nondeterministic finite state automaton</em> (NFA) from a regular
 *  expression and testing whether a given string is matched by that regular
 *  expression.
 *  It supports the following operations: <em>concatenation</em>,
 *  <em>closure</em>, <em>binary or</em>, and <em>parentheses</em>.
 *  It does not support <em>mutiway or</em>, <em>character classes</em>,
 *  <em>metacharacters</em> (either in the text or pattern),
 *  <em>capturing capabilities</em>, <em>greedy</em> or <em>relucantant</em>
 *  modifiers, and other features in industrial-strength implementations
 *  such as {@link java.util.regex.Pattern} and {@link java.util.regex.Matcher}.
 *  <p>
 *  This implementation builds the NFA using a digraph and a stack
 *  and simulates the NFA using digraph search (see the textbook for details).
 *  The constructor takes time proportional to <em>m</em>, where <em>m</em>
 *  is the number of characters in the regular expression.
 *  The <em>recognizes</em> method takes time proportional to <em>m n</em>,
 *  where <em>n</em> is the number of characters in the text.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/54regexp">Section 5.4</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class NFA {
  private Digraph nfa;
  private String pattern;
  public NFA(String pattern) {
    if (pattern == null) throw new IllegalArgumentException("null input");
    Deque<Integer> stack = new ArrayDeque<>();
    this.pattern = pattern;
    this.nfa = new Digraph(pattern.length() + 2);
    //the directed graph only records epsilon transitions

    for (int i = 1; i <= pattern.length(); i++) {
      char currentChar = pattern.charAt(i - 1);
      /* illegal syntax
       */
      if (currentChar == '*') {
        nfa.addEdge(i - 1, i);
        if (i - 2 >= 0 && pattern.charAt(i - 2) != ')') {
          nfa.addEdge(i, i - 1);
        }
        nfa.addEdge(i, i + 1);
      } else if (currentChar == '(') {
        stack.addFirst(i);
        nfa.addEdge(i - 1, i);
        nfa.addEdge(i, i + 1);
      } else if (currentChar == ')') {
        int lastAlternationIndex = -1;
        while(pattern.charAt(stack.peekFirst() - 1) != '(') {
          if (pattern.charAt(stack.peekFirst()-1) == '|') {
            int indexOnStack = stack.removeFirst();
            nfa.addEdge(indexOnStack - 1, i);
            if (lastAlternationIndex != -1) {
              nfa.addEdge(indexOnStack, lastAlternationIndex);
            }
            lastAlternationIndex = indexOnStack;
          }
        }
        if (pattern.charAt(stack.peekFirst()-1) != '(')
          throw new IllegalArgumentException("() must be paired");
        int leftParenthesisIndex = stack.removeFirst();
        if (i < pattern.length() && pattern.charAt(i) == '*') {
          nfa.addEdge(i + 1, leftParenthesisIndex);
        }
        if (lastAlternationIndex != -1) {
          nfa.addEdge(leftParenthesisIndex, lastAlternationIndex);
        }
        nfa.addEdge(i - 1, i);
        nfa.addEdge(i, i + 1);
      } else if (currentChar == '|') {
        stack.addFirst(i);
        nfa.addEdge(i, i + 1);
      } else {
        //regular chars or .
      }
    }
    //make sure start is connected to non-| chars
    if (pattern.length() >= 1 && pattern.charAt(0) != '|')
      nfa.addEdge(0, 1);
    if (!stack.isEmpty()) {
      int lastAlternationIndex = -1;
      while(!stack.isEmpty()) {
        int indexOnStack = stack.removeFirst();
        char charOnStack = pattern.charAt(indexOnStack - 1);
        if (charOnStack != '|')
          throw new IllegalArgumentException("non | found at the end of pattern, () not paired?");
        nfa.addEdge(indexOnStack - 1, nfa.V() - 1);
        if (lastAlternationIndex != -1) {
          nfa.addEdge(indexOnStack, lastAlternationIndex);
        }
        lastAlternationIndex = indexOnStack;
      }
      nfa.addEdge(0, lastAlternationIndex);
    }
  }
  public boolean recognizes(String q) {
    List<Integer> toVisit = new ArrayList<>();
    List<Integer> next = new ArrayList<>();

    DirectedDFS dfs = new DirectedDFS(nfa, 0);
    for (int j = 1; j < nfa.V(); j++) {
      if (dfs.marked(j)) {
        toVisit.add(j);
      }
    }
    if (dfs.marked(nfa.V() - 1))
      return true;

    for (int i = 0; i < q.length(); i++) {
      char currentChar = q.charAt(i);
      if (currentChar == '|' || currentChar == '*' || currentChar == '(' || currentChar == ')')
        throw new IllegalArgumentException("no | * ( ) allowed");
      for (int j : toVisit) {
        if (pattern.charAt(j -1) == '.' || pattern.charAt(j-1) == currentChar) {
          //if non-epsilon transition exists
          //2 nodes must be adjacent
          next.add(j + 1);
        }
      }
      dfs = new DirectedDFS(nfa, next);
      toVisit = new ArrayList<>();
      for (int j = 0; j < nfa.V(); j++) {
        if (dfs.marked(j)) {
          toVisit.add(j);
        }
      }
      if (dfs.marked(nfa.V() - 1))
        return true;
    }
    return false;
  }
  /**
   * Unit tests the {@code NFA} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    String regexp = "(" + args[0] + ")";
    String txt = args[1];
    NFA nfa = new NFA(regexp);
    StdOut.println(nfa.recognizes(txt));
  }
}

