package yunfeiImplementAlgs4;

import com.sun.org.apache.xpath.internal.operations.Div;

import edu.princeton.cs.algs4.Graph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

/******************************************************************************
 *  Compilation:  javac GraphClient.java
 *  Execution:    java GraphClient graph.txt
 *  Dependencies: Graph.java
 *
 *  Typical graph-processing code.
 *
 *  % java GraphClient tinyG.txt
 *  13 13
 *  0: 6 2 1 5 
 *  1: 0 
 *  2: 0 
 *  3: 5 4 
 *  4: 5 6 3 
 *  5: 3 4 0 
 *  6: 0 4 
 *  7: 8 
 *  8: 7 
 *  9: 11 10 12 
 *  10: 9 
 *  11: 9 12 
 *  12: 11 9 
 *
 *  vertex of maximum degree = 4
 *  average degree           = 2
 *  number of self loops     = 0
 *  
 ******************************************************************************/

public class GraphClient {
    //maximum degree
    public static int maxDegree(Graph G) {
        if (G == null) {
            throw new NullPointerException();
        }
        int maxDegree = -1;
        for (int i = 0; i < G.V(); i++) {
            maxDegree = Math.max(G.degree(i), maxDegree);            
        }
        return maxDegree;
    }
    //average degree
    public static int avgDegree(Graph G) {
        if (G == null) {
            throw new NullPointerException();
        }
        if (G.V() == 0) {
            throw new IllegalArgumentException("zero vertices");
        }
        //return G.E()/G.V(); //wrong!!!
        return 2*G.E()/G.V(); //each edge is incident on two vertices
    }
    //number of self-loops
    public static int numberOfSelfLoops(Graph G) {
        if (G == null) {
            throw new NullPointerException();
        }
        int selfLoopCount = 0;
        for (int i = 0; i < G.V(); i++) {
            for (int j : G.adj(i)) {
                if (j == i) {
                    selfLoopCount++;
                    break; //because we break out of the loop once we determine this is a self-loop
                            //there is no need to divide it by 2
                }
            }
        }  
        return selfLoopCount;
    }
    // number of self-loops version 2
    public static int numberOfSelfLoops2(Graph G) {
        int count = 0;
        for (int v = 0; v < G.V(); v++)
            for (int w : G.adj(v))
                if (v == w) count++;
        return count/2;   // self loop appears in adjacency list twice, why? check out Graph API addEdge
    } 
    public static void main(String[] args) {
        In in = new In(args[0]);
        Graph G = new Graph(in);
        StdOut.println(G);


        StdOut.println("vertex of maximum degree = " + maxDegree(G));
        StdOut.println("average degree           = " + avgDegree(G));
        StdOut.println("number of self loops     = " + numberOfSelfLoops(G));
        StdOut.println("number of self loops     = " + numberOfSelfLoops2(G));

    }
}
