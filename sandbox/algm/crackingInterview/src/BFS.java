import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.Hashtable;
import java.util.ArrayList;
import java.util.Queue;
import java.util.Deque;
public class BFS {
	/**
	 * BFS traversal
	 * @param root
	 * @return
	 */
	public static List<Integer> bfs(GraphNode root) {
		List<Integer> ret = new LinkedList<Integer>();
		Map<GraphNode, Boolean> visited = new Hashtable<GraphNode, Boolean>();
		Queue<GraphNode> toVisit = new LinkedList<GraphNode>();
		if(root == null) return(ret);
		toVisit.offer(root);
		while (! toVisit.isEmpty()) {
			int n = toVisit.size();
			for (int i = 0; i < n; i++) {
				GraphNode cur = toVisit.poll();			
				if(visited.containsKey(cur) && visited.get(cur)) {
					continue; //there could exist a cycle
				} else {
					for(GraphNode j:cur.neighbor) {
						if(j != null) toVisit.offer(j);
					}
					visited.put(cur, true);
					ret.add(cur.val);
				}
			}
		}
		return(ret);
	}
	public List<List<Integer>> levelOrder(TreeNode root) {
		List<List<Integer>> ret = new LinkedList<List<Integer>>();
		Queue<TreeNode> toVisit = new LinkedList<TreeNode>();
		if(root == null) return(ret);
		toVisit.offer(root);
		while (! toVisit.isEmpty()) {
			int n = toVisit.size();
			List<Integer> curL = new LinkedList<Integer>();
			for (int i = 0; i < n; i++) {
				TreeNode cur = toVisit.poll();			
				if (cur.left != null) toVisit.offer(cur.left);
				if (cur.right != null) toVisit.offer(cur.right);
				curL.add(cur.val);
			}
			ret.add(curL);
		}
		return(ret);
	}
	public static void main(String[] args) {
		GraphNode root = new GraphNode(0);
		GraphNode child1 = new GraphNode(1);
		GraphNode child2 = new GraphNode(2);
		child2.neighbor.add(root);
		root.neighbor.add(child1);
		root.neighbor.add(child2);
		child1.neighbor.add(new GraphNode(3));
		List<Integer> l = BFS.bfs(root);
		for(int i:l) {
			System.out.println(i);
		}
	}

}
