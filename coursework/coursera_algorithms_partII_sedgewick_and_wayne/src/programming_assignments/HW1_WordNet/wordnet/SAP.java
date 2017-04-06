package programming_assignments.HW1_WordNet.wordnet;
import edu.princeton.cs.algs4.BreadthFirstDirectedPaths;
import edu.princeton.cs.algs4.Digraph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;
public class SAP {
    private final Digraph G;
    // constructor takes a digraph (not necessarily a DAG)
    public SAP(Digraph G) {
        if (G == null) {
            throw new java.lang.NullPointerException();
        }
        this.G = clone(G); //we want SAP immutable, however, Digraph is mutable, so we need to create a copy of it       
    }
    /**
     * return a deep copy of Digraph G
     * @param G
     * @return
     */
    private Digraph clone(Digraph G) {
        Digraph copy = new Digraph(G.V());
        for (int i = 0; i < G.V(); i++) {
            for (int j : G.adj(i)) {
                copy.addEdge(i, j);
            }
        }
        return copy;
    }

    // length of shortest ancestral path between v and w; -1 if no such path
    public int length(int v, int w) {
        if ( v < 0 || v >= G.V() || w < 0 || w >= G.V()) {
            throw new java.lang.IndexOutOfBoundsException();
        }
        BreadthFirstDirectedPaths bfsV = new BreadthFirstDirectedPaths(G, v);
        BreadthFirstDirectedPaths bfsW = new BreadthFirstDirectedPaths(G, w);
        int ancestor = -1;
        int minLength = -1;
        for (int i = 0; i < G.V(); i++) {
            if (bfsV.hasPathTo(i) && bfsW.hasPathTo(i)) {
                if (ancestor == -1 || bfsV.distTo(i) + bfsW.distTo(i) < minLength) {
                    ancestor = i;
                    minLength = bfsV.distTo(i) + bfsW.distTo(i);
                }
            }
        }
        return minLength;
    }

    // a common ancestor of v and w that participates in a shortest ancestral path; -1 if no such path
    public int ancestor(int v, int w) {
        if ( v < 0 || v >= G.V() || w < 0 || w >= G.V()) {
            throw new java.lang.IndexOutOfBoundsException();
        }
        BreadthFirstDirectedPaths bfsV = new BreadthFirstDirectedPaths(G, v);
        BreadthFirstDirectedPaths bfsW = new BreadthFirstDirectedPaths(G, w);
        int ancestor = -1;
        int minLength = -1;
        for (int i = 0; i < G.V(); i++) {
            if (bfsV.hasPathTo(i) && bfsW.hasPathTo(i)) {
                if (ancestor == -1 || bfsV.distTo(i) + bfsW.distTo(i) < minLength) {
                    ancestor = i;
                    minLength = bfsV.distTo(i) + bfsW.distTo(i);
                }
            }
        }
        return ancestor;
    }    

    // length of shortest ancestral path between any vertex in v and any vertex in w; -1 if no such path
    public int length(Iterable<Integer> v, Iterable<Integer> w) {
        if (v == null || w == null) {
            throw new java.lang.NullPointerException();
        }
        BreadthFirstDirectedPaths bfsV = new BreadthFirstDirectedPaths(G, v);
        BreadthFirstDirectedPaths bfsW = new BreadthFirstDirectedPaths(G, w);
        int ancestor = -1;
        int minLength = -1;
        for (int i = 0; i < G.V(); i++) {
            if (bfsV.hasPathTo(i) && bfsW.hasPathTo(i)) {
                if (ancestor == -1 || bfsV.distTo(i) + bfsW.distTo(i) < minLength) {
                    ancestor = i;
                    minLength = bfsV.distTo(i) + bfsW.distTo(i);
                }
            }
        }
        return minLength;
    }

    // a common ancestor that participates in shortest ancestral path; -1 if no such path
    public int ancestor(Iterable<Integer> v, Iterable<Integer> w) {
        if (v == null || w == null) {
            throw new java.lang.NullPointerException();
        }
        BreadthFirstDirectedPaths bfsV = new BreadthFirstDirectedPaths(G, v);
        BreadthFirstDirectedPaths bfsW = new BreadthFirstDirectedPaths(G, w);
        int ancestor = -1;
        int minLength = -1;
        for (int i = 0; i < G.V(); i++) {
            if (bfsV.hasPathTo(i) && bfsW.hasPathTo(i)) {
                if (ancestor == -1 || bfsV.distTo(i) + bfsW.distTo(i) < minLength) {
                    ancestor = i;
                    minLength = bfsV.distTo(i) + bfsW.distTo(i);
                }
            }
        }
        return ancestor;
    }

    // do unit testing of this class
    public static void main(String[] args) {
        In in = new In(args[0]);
        Digraph G = new Digraph(in);
        SAP sap = new SAP(G);
        while (!StdIn.isEmpty()) {
            int v = StdIn.readInt();
            int w = StdIn.readInt();
            int length   = sap.length(v, w);
            int ancestor = sap.ancestor(v, w);
            StdOut.printf("length = %d, ancestor = %d\n", length, ancestor);
        }
    }
}
