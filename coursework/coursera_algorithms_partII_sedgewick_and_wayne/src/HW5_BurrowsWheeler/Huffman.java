package HW5_BurrowsWheeler;


import java.util.Queue;
import java.nio.channels.IllegalSelectorException;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;

import edu.princeton.cs.algs4.BinaryIn;
import edu.princeton.cs.algs4.BinaryOut;

/******************************************************************************
 *  Compilation:  javac Huffman.java
 *  Execution:    java Huffman - < input.txt   (compress)
 *  Execution:    java Huffman + < input.txt   (expand)
 *  Dependencies: BinaryIn.java BinaryOut.java
 *  Data files:   http://algs4.cs.princeton.edu/55compression/abra.txt
 *                http://algs4.cs.princeton.edu/55compression/tinytinyTale.txt
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

/**
 *  The <tt>Huffman</tt> class provides static methods for compressing
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
    /* roadmap
     * 
     * node for a trie
     * build trie
     * write trie
     * compress
     * read trie from bit stream
     * expand
     */
    private Huffman() {} //this class should not be instantiated
    private static class Node implements Comparable<Node> {
        /*frequency of the character or frequency
         * of all characters that are descendants of this node
         * assume it will not exceed 2^31-1
         * in production code, we may have to use long instead
         */
        private int frequency;
        private char c;
        private Node left;
        private Node right;

        public Node(char c, Node left, Node right, int frequency) {
            this.c = c;
            this.left = left;
            this.right = right;
            this.frequency = frequency;
        }
        public Node(char c) {
            this.c = c;
            this.left = null;
            this.right = null;
            this.frequency = -1;
        }
        public Node() {
            this.c = '\0';
            this.left = null;
            this.right = null;
            this.frequency = -1;
        }
        public boolean isLeafNode() {
            return this.left == null && this.right == null;
        }

        public int compareTo(Node that) {            
            return this.frequency - that.frequency;
        }
    }
    /**
     * read from stdin, build a trie, then return the root 
     * @return
     */
    private static Node buildTrie(BinaryIn binaryIn) {
        /* Huffman coding is a two-pass algorithm
         * in first pass, we count frequency of each character
         * construct a trie with least frequent character at the bottom
         * in 2nd pass, we convert each character to Huffman code
         */
        Map<Character, Integer> count = new HashMap<Character, Integer>();
        //StdIn has to be rest in 2nd pass
        while(!binaryIn.isEmpty()) {
            char c = binaryIn.readChar();
            if (count.containsKey(c)) {
                count.put(c, count.get(c) + 1);
            } else {
                count.put(c, 1);
            }
        }
        Queue<Node> q = new PriorityQueue<Node>();
        for (char c : count.keySet()) {
            Node n = new Node(c, null, null, count.get(c));
            q.offer(n);
        }
        //if there is only one character
        if (q.size() == 1) {
            if (count.containsKey('\0')) {
                q.offer(new Node('\0',null,null,0));
            } else {
                q.offer(new Node('\1',null,null,0));
            }
        }

        while(q.size() > 1) {
            Node left = q.poll();
            Node right = q.poll();
            Node parent = new Node('\0', left, right, left.frequency + right.frequency);
            q.offer(parent);
        }
        return q.poll();
    }
    /**
     * write the trie in binary format
     * 01A01B1C
     * write 0 for non-leaf node, 1 + the ASCII code otherwise
     * @param root
     */
    private static void writeTrie(Node root, BinaryOut out) {
        if (root.isLeafNode()) {
            out.write(true);
            out.write(root.c, 8);
            //System.out.println('1');
           //System.out.println(root.c);
        } else {
            out.write(false);
            //System.out.println('0');
            writeTrie(root.left, out);
            writeTrie(root.right, out);
        }
    }

    public static void compress(String file, String output) {
        BinaryOut out = new BinaryOut(output);
        Node root = buildTrie(new BinaryIn(file));
        writeTrie(root, out);
        out.write(root.frequency);
        BinaryIn in = new BinaryIn(file);                
        Map<Character, String> huffmanCode = new HashMap<Character, String>(); 
        getHuffmanCode(root, huffmanCode, new StringBuilder());
        while(!in.isEmpty()) {
            char c = in.readChar();
            String code = huffmanCode.get(c);
            for (int i = 0; i < code.length(); i++) {
                if (code.charAt(i) == '0') {
                    out.write(false);
                    //System.out.println('0');
                } else if (code.charAt(i) == '1') {
                    out.write(true);
                    //System.out.println('1');
                } else 
                    throw new IllegalSelectorException();
            }
        }
        out.close();
    }
    /** given Trie root
     * return mapping of character to binary string (implicitly)
     * @param root
     * @return
     */
    private static void getHuffmanCode(Node root, Map<Character, String> code, StringBuilder sb) {
        if (root.isLeafNode()) {
            code.put(root.c, sb.toString());
        } else {
            sb.append('0');//go left
            getHuffmanCode(root.left, code, sb);
            sb.deleteCharAt(sb.length() - 1);
            sb.append('1');
            getHuffmanCode(root.right, code, sb);
            sb.deleteCharAt(sb.length() - 1);
        }
    }
    /**
     * read bit stream, convert huffman code to unicode character
     * @param file
     */
    public static void expand(String file, String output) {
        BinaryIn in = new BinaryIn(file);
        BinaryOut out = new BinaryOut(output);
        Node root = readTrie(in);
        int length = in.readInt();
        for (int i = 0; i < length; i++) {
            Node currentNode = root;
            //every time, we begin with the root
            //go down one level at a time for every bit
            while(!currentNode.isLeafNode()) {                
                if (in.readBoolean())  {
                    currentNode = currentNode.right;
                } else {
                    currentNode = currentNode.left;
                }
            }
            out.write(currentNode.c);
        }   
        out.close();
    }
    private static Node readTrie(BinaryIn in) {        
        if (in.readBoolean()) {
            //this is a leaf node
            return new Node(in.readChar());
        } else {
            Node root = new Node();
            root.left = readTrie(in);
            root.right = readTrie(in);
            return root;
        }
    }
    public static void main(String[] args) {
        if (args[0].equals("-")) {
            compress(args[1], args[2]);
        } else if (args[0].equals("+")) {
            expand(args[1], args[2]);
        } else
            throw new UnsupportedOperationException();
    }
}
