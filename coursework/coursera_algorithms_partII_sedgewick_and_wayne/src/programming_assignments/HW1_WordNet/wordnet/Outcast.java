package programming_assignments.HW1_WordNet.wordnet;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

public class Outcast {
    private WordNet net;
   public Outcast(WordNet wordnet) {
       // constructor takes a WordNet object
       net = wordnet;
   }
   public String outcast(String[] nouns) {
       // given an array of WordNet nouns, return an outcast
       String result = null;
       int[][] lengthMatrix = new int[nouns.length][nouns.length];
       int maxDistance = 0;
       for (int i = 0; i < nouns.length; i++) {
           for (int j = 0; j < i; j++) {
               lengthMatrix[j][i] = lengthMatrix[i][j] = net.distance(nouns[i], nouns[j]);
           }
       }
       for (int i = 0; i < nouns.length; i++) {
           String current = nouns[i];
           int currentDistance = 0;
           for (int j = 0; j < nouns.length; j++) {                                  
               currentDistance += lengthMatrix[i][j];
           }
           if (currentDistance > maxDistance) {
               result = current;
               maxDistance = currentDistance;
           }
       }
       return result;
   }
   public static void main(String[] args) {
       // see test client below
       WordNet wordnet = new WordNet(args[0], args[1]);
       Outcast outcast = new Outcast(wordnet);
       for (int t = 2; t < args.length; t++) {
           In in = new In(args[t]);
           String[] nouns = in.readAllStrings();
           StdOut.println(args[t] + ": " + outcast.outcast(nouns));
       }
   }
}
