package yunfeiImplementAlgs4;

import java.util.ArrayDeque;
import java.util.Deque;

import edu.princeton.cs.algs4.Digraph;

public class NFA {
    private Digraph eTransition; //episilon transition
    private char[] regex;
    private int regexLength;
    public NFA(String regexStr) {
        this.regex = regexStr.toCharArray();
        this.regexLength = regex.length;
        
        Deque<Integer> deque = new ArrayDeque<Integer>();
        eTransition = new Digraph(regexLength + 1); 
        for (int i = 0; i < regexLength; i++) {
            int lp = i;
            if (regex[i] == '(' || regex[i] == '|') {
                deque.addLast(i);
            } else if (regex[i] == ')') {
                int or = deque.removeLast();
                if (regex[or] == '|') {
                    lp = deque.removeLast();
                    eTransition.addEdge(or, i);
                    eTransition.addEdge(lp, or + 1);                   
                } else
                    lp = or;                
            }
            //look ahead for * closure
            //closure only appears after right parenthesis or 
            //a regular character
            if (i < regexLength - 1 && regex[i + 1] == '*') {
                eTransition.addEdge(lp, i + 1);
                eTransition.addEdge(i + 1, lp);
            }
            //add e-transition for (,),*, but not for or operator
            //when i == regexLength - 1, add one additional
            //e-transition to accept state
            if (regex[i] == '(' || regex[i] == ')' || regex[i] == '*') {
                eTransition.addEdge(i, i + 1);
            }
        }
    }
    public boolean recognizes(String txt) {
        Bag<Integer> next = new Bag<Integer>();
        DirectedDFS dfs = new DirectedDFS(eTransition, 0);
        
        //first add all nodes directly reachable
        //to search space
        for (int i = 0; i < eTransition.V(); i++) {
            if (dfs.marked(i)) {
                next.add(i);
                if (i == regexLength) {
                    return true;
                }
            }            
        }
        for (int i = 0; i < txt.length(); i++) {
            char current = txt.charAt(i);
            Bag<Integer> matches = new Bag<Integer>();
            for (int j : next) {                
                if (regex[j] == '.' || regex[j] == current) {
                    matches.add(j + 1); //since we only store epsilon transitions, we must
                                        //add next adjacent char index to search space
                }
            }
            next = new Bag<Integer>();
            dfs = new DirectedDFS(eTransition, matches);
            for (int k = 0; k < eTransition.V(); k++) {                
                if (dfs.marked(k)) {
                    next.add(k);
                    if (k == regexLength) {
                        return true;
                    }
                }
            }                        
        }
        return false;
    }
    public static void main(String[] args) {
        NFA test = new NFA("AB((C|D*E)F)*G");
        System.out.println(test.recognizes("ABCFDDEFG"));
        System.out.println(test.recognizes("ABCCFDDEFG") == false);
    }
}
