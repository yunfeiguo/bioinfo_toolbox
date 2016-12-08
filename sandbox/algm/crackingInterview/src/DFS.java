import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedList;
public class DFS {
	public static List<Integer> dfs(GraphNode root) {
		List<Integer> ret = new LinkedList<Integer>();
		Map<GraphNode, Boolean> visited = new Hashtable<GraphNode, Boolean>();
		Stack<GraphNode> toVisit = new Stack<GraphNode>();
		if(root == null) return(ret);
		toVisit.push(root);
		while (!toVisit.isEmpty()) {			
			GraphNode cur = toVisit.pop();
			if (visited.containsKey(cur) && visited.get(cur)) {
				continue;
			} else {
				visited.put(cur, true);
				for (GraphNode j:cur.neighbor) {
					toVisit.push(j);
				}
				ret.add(cur.val);
			}
		}
		Collections.reverse(ret);
		return(ret);
	}
	public static List<Integer> dfsRecur(GraphNode root, Map<GraphNode, Boolean> visited) {
		List<Integer> ret = new LinkedList<Integer>();
		if (root == null) return(ret);
		visited.put(root, true);
		if (root.neighbor != null && !root.neighbor.isEmpty()) {
			for(GraphNode i:root.neighbor) {
				if (visited.containsKey(i) && visited.get(i)) {
					continue;
				} else {
					ret.addAll(DFS.dfsRecur(i, visited));
				}
			}
		}
		ret.add(root.val);
		return(ret);
	}
	public static void main(String[] args) {
		Map<GraphNode, Boolean> visited = new Hashtable<GraphNode, Boolean>();
		GraphNode root = new GraphNode(0);
		GraphNode child1 = new GraphNode(1);
		GraphNode child2 = new GraphNode(2);
		child2.neighbor.add(root);
		root.neighbor.add(child1);
		root.neighbor.add(child2);
		child1.neighbor.add(new GraphNode(3));
		//List<Integer> l = DFS.dfs(root);
		List<Integer> l = DFS.dfsRecur(root, visited);
		for(int i:l) {
			System.out.println(i);
		}
	}
}
