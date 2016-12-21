package yunfeiImplementAlgs4;
/******************************************************************************
 *  Compilation:  javac SymbolGraph.java
 *  Execution:    java SymbolGraph filename.txt delimiter
 *  Dependencies: ST.java Graph.java In.java StdIn.java StdOut.java
 *  Data files:   http://algs4.cs.princeton.edu/41graph/routes.txt
 *                http://algs4.cs.princeton.edu/41graph/movies.txt
 *                http://algs4.cs.princeton.edu/41graph/moviestiny.txt
 *                http://algs4.cs.princeton.edu/41graph/moviesG.txt
 *                http://algs4.cs.princeton.edu/41graph/moviestopGrossing.txt
 *  
 *  %  java SymbolGraph routes.txt " "
 *  JFK
 *     MCO
 *     ATL
 *     ORD
 *  LAX
 *     PHX
 *     LAS
 *
 *  % java SymbolGraph movies.txt "/"
 *  Tin Men (1987)
 *     Hershey, Barbara
 *     Geppi, Cindy
 *     Jones, Kathy (II)
 *     Herr, Marcia
 *     ...
 *     Blumenfeld, Alan
 *     DeBoy, David
 *  Bacon, Kevin
 *     Woodsman, The (2004)
 *     Wild Things (1998)
 *     Where the Truth Lies (2005)
 *     Tremors (1990)
 *     ...
 *     Apollo 13 (1995)
 *     Animal House (1978)
 *
 * 
 *  Assumes that input file is encoded using UTF-8.
 *  % iconv -f ISO-8859-1 -t UTF-8 movies-iso8859.txt > movies.txt
 *
 ******************************************************************************/

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/**
 *  The <tt>SymbolGraph</tt> class represents an undirected graph, where the
 *  vertex names are arbitrary strings.
 *  By providing mappings between string vertex names and integers,
 *  it serves as a wrapper around the
 *  {@link Graph} data type, which assumes the vertex names are integers
 *  between 0 and <em>V</em> - 1.
 *  It also supports initializing a symbol graph from a file.
 *  <p>
 *  This implementation uses an {@link ST} to map from strings to integers,
 *  an array to map from integers to strings, and a {@link Graph} to store
 *  the underlying graph.
 *  The <em>index</em> and <em>contains</em> operations take time 
 *  proportional to log <em>V</em>, where <em>V</em> is the number of vertices.
 *  The <em>name</em> operation takes constant time.
 *  <p>
 *  For additional documentation, see <a href="http://algs4.cs.princeton.edu/41graph">Section 4.1</a> of
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class SymbolGraph {
    private Graph G; //the Graph obj
    private int V; //# of vertices
    private Map<String, Integer> name2index; //store name->index mapping
    private ArrayList<String> index2name; //store index->name mapping    
    /**  
     * Initializes a graph from a file using the specified delimiter.
     * Each line in the file contains
     * the name of a vertex, followed by a list of the names
     * of the vertices adjacent to that vertex, separated by the delimiter.
     * @param filename the name of the file
     * @param delimiter the delimiter between fields
     */
    public SymbolGraph(String filename, String delimiter) {
        name2index = new HashMap<String, Integer>();
        index2name = new ArrayList<String>();
        V = 0;
        LineReader r = new LineReader(filename);
        for (String s : r) {
            addVertices(s.split(delimiter));
        }
        r = new LineReader(filename);
      //second time, we add all edges
        G = new Graph(V);
        for (String s : r) {
            addEdges(s.split(delimiter));
        }             
    }    
    private void addVertices(String[] a) {
        if (a == null || a.length < 1) {
            return;
        }
        for (int i = 0; i < a.length; i++) {
            if (!contains(a[i])) {
                name2index.put(a[i], V++);                
                index2name.add(a[i]);
            }
        }        
    }
    private void addEdges(String[] a) {
        if (a == null || a.length < 1) {
            return;
        }
        String vName = a[0];
        int v = index(vName);
        for (int i = 1; i < a.length; i++) {
            G.addEdge(v, index(a[i]));
        }
    }
    /**
     * Does the graph contain the vertex named <tt>s</tt>?
     * @param s the name of a vertex
     * @return <tt>true</tt> if <tt>s</tt> is the name of a vertex, and <tt>false</tt> otherwise
     */
    public boolean contains(String s) {
        return name2index.containsKey(s);
    }
    /**
     * Returns the integer associated with the vertex named <tt>s</tt>.
     * @param s the name of a vertex
     * @return the integer (between 0 and <em>V</em> - 1) associated with the vertex named <tt>s</tt>
     */
    public int index(String s) {
        return name2index.get(s);
    }
    /**
     * Returns the name of the vertex associated with the integer <tt>v</tt>.
     * @param v the integer corresponding to a vertex (between 0 and <em>V</em> - 1) 
     * @return the name of the vertex associated with the integer <tt>v</tt>
     */
    public String name(int v) {
        return index2name.get(v);
    }
    
    /**
     * Returns the graph assoicated with the symbol graph. It is the client's responsibility
     * not to mutate the graph.
     * @return the graph associated with the symbol graph
     */
    public Graph G() {
        return G;
    }
    /**
     * Unit tests the <tt>SymbolGraph</tt> data type.
     */
    public static void main(String[] args) {
        String filename  = args[0];
        String delimiter = args[1];
        SymbolGraph sg = new SymbolGraph(filename, delimiter);
        Graph G = sg.G();
        while (StdIn.hasNextLine()) {
            String source = StdIn.readLine();
            if (sg.contains(source)) {
                int s = sg.index(source);
                for (int v : G.adj(s)) {
                    StdOut.println("   " + sg.name(v));
                }
            }
            else {
                StdOut.println("input not contain '" + source + "'");
            }
        }
    }
}
