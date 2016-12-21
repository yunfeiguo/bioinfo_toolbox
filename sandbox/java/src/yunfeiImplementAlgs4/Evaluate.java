/******************************************************************************
 *  Compilation:  javac Evaluate.java
 *  Execution:    java Evaluate
 *  Dependencies: Stack.java
 *
 *  Evaluates (fully parenthesized) arithmetic expressions using
 *  Dijkstra's two-stack algorithm.
 *
 *  % java Evaluate 
 *  ( 1 + ( ( 2 + 3 ) * ( 4 * 5 ) ) ) 
 *  101.0 
 *
 *  % java Evaulate
 *  ( ( 1 + sqrt ( 5 ) ) / 2.0 ) 
 *  1.618033988749895
 *
 *
 *  Note: the operators, operands, and parentheses must be
 *  separated by whitespace. Also, each operation must
 *  be enclosed in parentheses. For example, you must write
 *  ( 1 + ( 2 + 3 ) ) instead of ( 1 + 2 + 3 ).
 *  See EvaluateDeluxe.java for a fancier version.
 *
 *
 *  Remarkably, Dijkstra's algorithm computes the same
 *  answer if we put each operator *after* its two operands
 *  instead of *between* them.
 *
 *  % java Evaluate
 *  ( 1 ( ( 2 3 + ) ( 4 5 * ) * ) + ) 
 *  101.0
 *
 *  Moreover, in such expressions, all parentheses are redundant!
 *  Removing them yields an expression known as a postfix expression.
 *  1 2 3 + 4 5 * * + 
 * 
 *
 ******************************************************************************/
package yunfeiImplementAlgs4;

import java.util.Deque;
import java.util.LinkedList;

import edu.princeton.cs.algs4.StdIn;

public class Evaluate {
    public static void main(String[] args) {
        Deque<String> op = new LinkedList<String>();
        Deque<Double> val = new LinkedList<Double>();

        while (!StdIn.isEmpty()) {
            String s = StdIn.readString();
            if (s.equals("(")) ;
            else if (s.equals("+")) op.addFirst("+");
            else if (s.equals("*")) op.addFirst("*");
            else if (s.equals("-")) op.addFirst("-");
            else if (s.equals("/")) op.addFirst("/");
            else if (s.equals(")")) {
                String operator = op.removeFirst();
                Double n1 = val.removeFirst();
                Double n2 = val.removeFirst();
                if (operator.equals("+")) val.addFirst(n1 + n2);
                else if (operator.equals("*")) val.addFirst(n1 * n2);
                else if (operator.equals("-")) val.addFirst(n2 - n1);
                else if (operator.equals("/")) val.addFirst(n2 / n1);
            }
            else val.addFirst(Double.parseDouble(s));      
        }
        System.out.println(val.removeFirst());
    }

}
