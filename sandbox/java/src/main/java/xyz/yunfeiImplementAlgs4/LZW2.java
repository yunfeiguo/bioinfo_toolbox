package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;
import edu.princeton.cs.algs4.TST;

/******************************************************************************
 *  Compilation:  javac LZW2.java
 *  Execution:    java LZW2 - input.txt output.txt.lzw  (compress)
 *  Execution:    java LZW2 + input.txt.lzw output.txt   (expand)
 *  Dependencies: BinaryIn.java BinaryOut.java
 *
 *  Compress or expand binary input from standard input using LZW2.
 *
 *  WARNING: STARTING WITH ORACLE JAVA 6, UPDATE 7 the SUBSTRING
 *  METHOD TAKES TIME AND SPACE LINEAR IN THE SIZE OF THE EXTRACTED
 *  SUBSTRING (INSTEAD OF CONSTANT SPACE AND TIME AS IN EARLIER
 *  IMPLEMENTATIONS).
 *
 *  See <a href = "http://java-performance.info/changes-to-string-java-1-7-0_06/">this article</a>
 *  for more details.
 *
 ******************************************************************************/

/**
 *  The <tt>LZW2</tt> class provides static methods for compressing
 *  and expanding a binary input using LZW2 compression over the 8-bit extended
 *  ASCII alphabet with 12-bit codewords.
 *  <p>
 *  For additional documentation,
 *  see <a href="http://algs4.cs.princeton.edu/55compress">Section 5.5</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick  
 *  @author Kevin Wayne
 */
public class LZW2 {
    private static int R = 256; //dictionary size for uncompressed characters
    private static int L = 4096; //dictionary size for compressed characters
    private static int lgL = 12; //bits of L
    public static void compress(String input, String output) {
        TST<Integer> trie = new TST<Integer>();
        //put single byte characters into trie
        int i;
        for (i = 0; i < R; i++) {
            trie.put(""+(char) i, i);
        }
        i = R + 1; //R is reserved for end of file
        BinaryIn in = new BinaryIn(input);
        BinaryOut out = new BinaryOut(output);
        StringBuilder sb = new StringBuilder();
        while(!in.isEmpty()) {
            Character c = in.readChar();
            sb.append(c);
            if (trie.contains(sb.toString())) {
                continue;
            } else {
                if (i < L) {
                    trie.put(sb.toString(), i++);
                }
                sb.deleteCharAt(sb.length() - 1);
                //there is a performance penalty to repeatedly call toString
                //and trie.get(), this can be improved
                //by implementing a trie inside this class
                //and go down one node at a time
                out.write(trie.get(sb.toString()), lgL);
                sb = new StringBuilder();
                sb.append(c);
            }            
        }
        if (sb.length() > 0) {
            out.write(trie.get(sb.toString()), lgL);
        }
        out.write(R, lgL); //write EOF
        out.close();
    }
    public static void expand(String input, String output) {
           BinaryIn in = new BinaryIn(input);
           BinaryOut out = new BinaryOut(output);
           String[] encoding = new String[L];
           int i = 0;
           //store all single-byte chars
           for (i = 0; i < R; i++) {
               encoding[i] = ""+(char) i;
           }
           i = R + 1; //skip R (EOF)
           String current = "";
           while(true) {
               int code = in.readInt(lgL);
               if (code == R) {//EOF
                   break;
               }
               if (code == i) {
                   //if the code is exactly next code we 
                   //want to add, then we know for sure
                   //the first char in previous string
                   //is same as the first char in current
                   //string
                   current += current.charAt(0);        
                   out.write(current);
               } else {
                   out.write(encoding[code]);
                   current += encoding[code].charAt(0);                   
               }
               if (i < L && current.length() > 1) {
                   //>1 is necessary otherwise we will waste
                   //space storing single-character strings
                   encoding[i++] = current;
               }
               current = encoding[code];
           }
           out.close();
    }
    public static void main(String[] args) {
        if (args[0].equals("-")) {
            compress(args[1], args[2]);
        } else if (args[0].equals("+")) {
            expand(args[1], args[2]);
        } else {
            throw new IllegalArgumentException();
        }
    }
}
