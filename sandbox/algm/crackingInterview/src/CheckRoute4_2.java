import java.util.*;
public class CheckRoute4_2 {
	/**
	 * check whether a route exists from s to e
	 * @param s
	 * @param e
	 * @return
	 */
	public static boolean check (GraphNode s, GraphNode e) {
		//use BFS
		Map<GraphNode, Boolean> visited = new Hashtable<GraphNode, Boolean>();
		Queue<GraphNode> toVisit = new LinkedList<GraphNode>();
		if(s == null || e == null) return(false);
		toVisit.offer(s);
		while(!toVisit.isEmpty()) {
			int n = toVisit.size();
			for (int i = 0; i < n; i++) {
				GraphNode v = toVisit.poll();
				if(visited.containsKey(v) && visited.get(v)) continue;
				if (v == e) return(true);
				toVisit.addAll(v.neighbor);
				visited.put(v, true);
			}
		}
		return(false);		
	}
	public static boolean checkDFS (GraphNode s, GraphNode e) {
		//use DFS		
		if(s == null || e == null) return(false);
		Map<GraphNode, Boolean> visited = new Hashtable<GraphNode, Boolean>();
		Stack<GraphNode> stack = new Stack<GraphNode>();
		stack.push(s);
		while(!stack.isEmpty()) {
			GraphNode cur = stack.pop();
			if(visited.containsKey(cur) && visited.get(cur)) continue;
			if (cur == e) return(true);
			for(GraphNode i:cur.neighbor) stack.push(i);
			visited.put(cur, true);
		}
		return(false);
	}
    public static List<String> binaryTreePaths(TreeNode root) {
        List<String> ret = new LinkedList<String>();
        if (root == null) return(ret);
        Stack<TreeNode> s = new Stack<TreeNode>();
        Stack<String> curPath = new Stack<String>(); //record path from root to current node
        Map<TreeNode, Boolean> visited = new Hashtable<TreeNode, Boolean>();
        s.push(root);
        curPath.push(""+root.val);
        while(!s.isEmpty()) {
            TreeNode cur = s.peek();
            if(cur.left == null && cur.right == null) {
                //a leaf node
                ret.add(concatStr(curPath));
                curPath.pop();//pop out the leaf node from current paths
                s.pop();//pop out the leaf from nodes to be visited
            } else {
            //not a leaf node
            if (cur.left != null) {
                if (!(visited.containsKey(cur.left) && visited.get(cur.left))) {
                                 curPath.push(""+cur.left);
                s.push(cur.left);
                continue;   
                }
            }
            if (cur.right != null) {
                if (!(visited.containsKey(cur.left) && visited.get(cur.left))) {
                curPath.push(""+cur.right);
                s.push(cur.right);
                continue;
                }
            }
            curPath.pop();
            s.pop();
            }
            visited.put(cur,true);//this node has been visited
        }
        return(ret);
    }
    private static String concatStr(Stack<String> s) {
        String ret = new String();
        int n = s.size();
        for (String i:s) {
        	if (n == 1) {
        		ret = ret + i;
        	} else {
        		ret = ret + i + "->";
        	}
        	n--;
        }
        return(ret);
    }
	public static void main(String[] args) {
		GraphNode s = new GraphNode(1);
		System.out.println(binaryTreePaths(new TreeNode(1)));
		GraphNode e = new GraphNode(2);
		e.neighbor.add(s);
		System.out.println("This should return false");
		System.out.println(CheckRoute4_2.check(s,e));
		System.out.println(CheckRoute4_2.checkDFS(s,e));
		s.neighbor.add(e);
		System.out.println("This should return true");
		System.out.println(CheckRoute4_2.check(s,e));
		System.out.println(CheckRoute4_2.checkDFS(s,e));
		Stack<String> x = new Stack<String>();
		x.push("1");
		x.push("2");
		System.out.println(concatStr(x));
	}

}
