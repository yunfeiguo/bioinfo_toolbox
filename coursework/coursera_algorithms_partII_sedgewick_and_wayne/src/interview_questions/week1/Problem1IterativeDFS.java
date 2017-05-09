package interview_questions.week1;
import edu.princeton.cs.algs4.DepthFirstSearch;
import edu.princeton.cs.algs4.Graph;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;
import org.apache.log4j.Logger;

import java.util.*;

/**
 problem:
 Nonrecursive depth-first search.
 Implement depth-first search in an undirected graph without using recursion.

 answer:
 find all nodes connected to a source vertex in a graph
 traversal must be using iterative DFS
 */
public class Problem1IterativeDFS {
	final static Logger logger = Logger.getLogger(Problem1IterativeDFS.class);
	private Graph g;
	private int source;
	private boolean[] visited;
	private int count;
	//we can also record paths, but here I just ignore them
	public Problem1IterativeDFS(Graph g, int v) {
		this.g = g; //to make the program more robust, we need to make a copy of graph as it is mutable.
		this.source = v;
		this.visited = new boolean[g.V()];
		count = 0;
		dfs();
	}
	public boolean isConnected(int v) {
		return visited[v];
	}
	public int count() {
		return count;
	}
	private void dfs() {
		Deque<Iterator<Integer>> neighborIteratorStack = new ArrayDeque<>(); //record which neighbor to visit next
		visited[source] = true;
		count++;
		neighborIteratorStack.addFirst(g.adj(source).iterator());
		while(!neighborIteratorStack.isEmpty()) {
			Iterator<Integer> neighborIterator = neighborIteratorStack.removeFirst();
			while (neighborIterator.hasNext()) {
				int next = neighborIterator.next();
				if (visited[next])
					continue;
        visited[next] = true;
				count++;
				neighborIteratorStack.addFirst(neighborIterator);
				neighborIteratorStack.addFirst(g.adj(next).iterator());
				break;
			}
		}
	}

	/**
	 * Unit tests the {@code DepthFirstSearch} data type.
	 *
	 * @param args the command-line arguments
	 */
	public static void main(String[] args) {
    //simple test
		Graph g = new Graph(5);
		g.addEdge(0,1);
		g.addEdge(0,3);
		g.addEdge(0,2);
		g.addEdge(2,3);
		g.addEdge(2,4);
		Problem1IterativeDFS test = new Problem1IterativeDFS(g, 0);
		System.out.println(test.isConnected(1));
		System.out.println(test.isConnected(3));

    //more complex test
		In in = new In(args[0]);
		Graph G = new Graph(in);
		int s = Integer.parseInt(args[1]);
		logger.info("recursive DFS");
		DepthFirstSearch search = new DepthFirstSearch(G, s);
		for (int v = 0; v < G.V(); v++) {
			if (search.marked(v))
				StdOut.print(v + " ");
		}

		StdOut.println();
		if (search.count() != G.V()) logger.info("NOT connected");
		else                         logger.info("connected");
    //my iterative DFS
		logger.info("iterative DFS");
		Problem1IterativeDFS search2 = new Problem1IterativeDFS(G, s);
		for (int v = 0; v < G.V(); v++) {
			if (search2.isConnected(v))
				StdOut.print(v + " ");
		}

		StdOut.println();
		if (search2.count() != G.V()) logger.info("NOT connected");
		else                         logger.info("connected");
	}
}
