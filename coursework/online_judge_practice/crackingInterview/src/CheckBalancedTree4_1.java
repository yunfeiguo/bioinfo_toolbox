import java.util.Queue;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
public class CheckBalancedTree4_1 {
	public static boolean checkRecur(TreeNode root) {
		//recursive solution
		return(maxDepth(root) - minDepth(root) <= 1);
	}
	private static int maxDepth(TreeNode root) {
		if (root == null) return(0);
		return(Math.max(maxDepth(root.left), maxDepth(root.right)) + 1);
	}
	private static int minDepth(TreeNode root) {
		if (root == null) return(0);
		return(Math.min(minDepth(root.left), minDepth(root.right)) + 1);
	}
	public static boolean check(TreeNode root) {
		//use BFS to traverse the tree
		//output levels for leaf nodes
		//check if levels of leaf nodes differ by more than one
		//we should pay attention to the definition of balance here
		//is it defined by distance of leaf node to root node or by depth of subtrees?
		int lvl = 0;
		Queue<TreeNode> toVisit = new LinkedList<TreeNode>();
		List<Integer> allLvl = new LinkedList<Integer>();
		if (root == null) return(true);
		toVisit.offer(root);
		while (!toVisit.isEmpty()) {
			int n = toVisit.size();
			for(int i = 0; i < n; i++) {
				TreeNode cur = toVisit.poll();
				if(cur.left == null || cur.right == null) {
					//we need to record its level
					allLvl.add(lvl);
				} 	
				if (cur.left != null) {
					toVisit.offer(cur.left);
				}
				if (cur.right != null) {
					toVisit.offer(cur.right);
				}
			}
			lvl++;
		}
		
		/*
		for (Integer i:allLvl) {
			System.out.println(i);
		}*/
		return(Collections.max(allLvl)-Collections.min(allLvl) <= 1);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(1);				
		System.out.println("This should return true");
		System.out.println(CheckBalancedTree4_1.check(root));
		System.out.println(CheckBalancedTree4_1.checkRecur(root));
		root.left = new TreeNode(2);
		root.left.left = new TreeNode(3);
		System.out.println("This should return false");
		System.out.println(CheckBalancedTree4_1.check(root));	
		System.out.println(CheckBalancedTree4_1.checkRecur(root));
	}

}
