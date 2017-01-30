package HW1_WordNet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.DirectedCycle;
import edu.princeton.cs.algs4.In;
public class WordNet {    
    private final Digraph G;
    private final Map<String, List<Integer>> dictionary;
    private final ArrayList<String> synsetList;
    private final SAP SAPFinder;

    // constructor takes the name of the two input files
    public WordNet(String synsets, String hypernyms) {
        if (synsets == null || hypernyms == null) {
            throw new java.lang.NullPointerException();
        }
        //process input
        dictionary = new HashMap<String, List<Integer>>();
        synsetList = new ArrayList<String>();
        processSynsets(synsets);
        G = new Digraph(synsetList.size());
        processHypernyms(hypernyms);

        //sanity check
        DirectedCycle cycleFinder = new DirectedCycle(G);
        if (cycleFinder.hasCycle()) {
            throw new java.lang.IllegalArgumentException("there is cycle");
        }        
        boolean foundOneRoot = false;
        for (int i = 0; i < G.V(); i++) {
            if (G.outdegree(i) == 0) {
                if (foundOneRoot) {
                    throw new java.lang.IllegalArgumentException("multiple roots");
                } else {
                    foundOneRoot = true;
                }
            }
        }
        SAPFinder = new SAP(G);
    }    
    private void processHypernyms(String hypernyms) {
        In input = new In(hypernyms);
        while (input.hasNextLine()) {
            String[] fields = input.readLine().split(",");
            int hyponym = Integer.parseInt(fields[0]);
            if (hyponym < 0 || hyponym >= G.V()) {
                throw new java.lang.IllegalArgumentException();
            }
            for (int i = 1; i < fields.length; i++) {
                int hypernym = Integer.parseInt(fields[i]);
                if (hypernym < 0 || hypernym >= G.V()) {
                    throw new java.lang.IllegalArgumentException();
                }
                G.addEdge(hyponym, hypernym);
            }
        }
    }
    private void processSynsets(String synsets) {
        In input = new In(synsets);                
        while(input.hasNextLine()) {
            //example input
            //0,a,one
            //assume index is increasing            
            String[] fields = input.readLine().split(",");
            if (Integer.parseInt(fields[0]) == synsetList.size()) {                
                String[] words = fields[1].split(" ");
                //what if a word appears on multiple lines?                
                for (int i = 0; i < words.length; i++) {
                    if (!dictionary.containsKey(words[i])) {
                        dictionary.put(words[i], new ArrayList<Integer>());
                    }
                    dictionary.get(words[i]).add(synsetList.size());
                }
                synsetList.add(fields[1]);
            } else {
                throw new java.lang.IllegalArgumentException("index not increasing from 0 to G.V - 1");
            }            
        }                
    }

    // returns all WordNet nouns
    public Iterable<String> nouns() {
        return dictionary.keySet();
    }

    // is the word a WordNet noun?
    public boolean isNoun(String word) {
        if (word == null) {
            throw new java.lang.NullPointerException();
        }
        return dictionary.containsKey(word);
    }

    // distance between nounA and nounB (defined below)
    public int distance(String nounA, String nounB) {
        if (nounA == null || nounB == null) {
            throw new java.lang.NullPointerException();
        }
        if (isNoun(nounA) == false || isNoun(nounB) == false) {
            throw new java.lang.IllegalArgumentException();
        }        

        return SAPFinder.length(dictionary.get(nounA), dictionary.get(nounB));
    }

    // a synset (second field of synsets.txt) that is the common ancestor of nounA and nounB
    // in a shortest ancestral path (defined below)
    public String sap(String nounA, String nounB) {
        if (nounA == null || nounB == null) {
            throw new java.lang.NullPointerException();
        }
        if (isNoun(nounA) == false || isNoun(nounB) == false) {
            throw new java.lang.IllegalArgumentException();
        }        
        int ancestorID = SAPFinder.ancestor(dictionary.get(nounA), dictionary.get(nounB));
        return synsetList.get(ancestorID);
    }

    // do unit testing of this class
    public static void main(String[] args) {
        WordNet test = new WordNet(args[0], args[1]);
        System.out.println(test.sap("b", "c"));

    }
}
