package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 4/14/17.
 */

import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac Validate.java
 *  Execution:    java Validate pattern text
 *  Dependencies: StdOut.java
 *
 *  Reads in a regular expression and a text input string from the
 *  command line and prints true or false depending on whether
 *  the pattern matches the text.
 *
 *  % java Validate "..oo..oo." bloodroot
 *  true
 *
 *  % java Validate "..oo..oo." nincompoophood
 *  false
 *
 *  % java Validate "[^aeiou]{6}" rhythm
 *  true
 *
 *  % java Validate "[^aeiou]{6}" rhythms
 *  false
 *
 ******************************************************************************/

public class Validate {

  public static void main(String[] args) {
    //String regexp = args[0];
    //String text   = args[1];
    String regexp = "(a|aa)*b";
    for (String text : new String[]{"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaac",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaac ",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaac",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaac",
            }) {

      StdOut.println("backtracking");
      StdOut.print(regexp + " matching " + text + " ");
      long start = System.nanoTime();
      StdOut.print(text.matches(regexp));
      long end = System.nanoTime();
      StdOut.println(" " + (end-start)/1e9 + " seconds.");

      StdOut.println("no backtracking");
      StdOut.print(regexp + " matching " + text + " ");
      start = System.nanoTime();
      NFA nfa = new NFA(regexp);
      StdOut.print(nfa.recognizes(text));
      end = System.nanoTime();
      StdOut.println(" " + (end-start)/1e9 + " seconds.");
    }
  }

}
