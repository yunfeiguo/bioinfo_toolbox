import java.util.*;
public class Subtree4_7 {
	public static void main(String[] args) {
		/*
		 * when dealing with this question, first we need to set up appropriate assumptions
		 * how to check if two nodes are equal? by reference, by val?
		 * are there duplicate nodes?
		 */
		TreeNode root = new TreeNode(0);
		Stack<TreeNode> s = new Stack<TreeNode>();
		TreeNode subtree = null;
		int count = 1000;
		s.push(root);
		while(s.size()>0 && count > 0) {
			TreeNode cur = s.pop();
			if (count == 100) {
				subtree = cur;
			}
			cur.left = new TreeNode((int) Math.random()*10000000);
			s.push(cur.left);
			cur.right = new TreeNode((int) Math.random()*10000000);
			s.push(cur.right);
			count--;
		}
		System.out.println(subtree(root, subtree));
		System.out.println(subtree2(root, subtree));
	}
	/**
	 * determine if subtree is a subtree of root
	 * @param root
	 * @param subtree
	 * @return
	 */
	public static boolean subtree(TreeNode root, TreeNode subtree) {
		boolean foundSubtreeRoot = false;
		Stack<TreeNode> s = new Stack<TreeNode>();
		Queue<TreeNode> subtreePreorder = preorderTraversal(subtree);
		if (root == null || subtree == null) return(false);
		s.push(root);
		while(s.size() > 0) {
			TreeNode cur = s.pop();
			if(!foundSubtreeRoot) {
				if (cur.val == subtree.val) {
					foundSubtreeRoot = true;
					subtreePreorder.poll();
				}
			} else {
				if(subtreePreorder.isEmpty()) return(true);
				TreeNode subtreeNext = subtreePreorder.poll();
				if(subtreeNext.val != cur.val) return(false);
			}
			if(cur.right!=null) s.push(cur.right);
			if(cur.left!=null) s.push(cur.left);
		}
		return(false);
	}
	public static Queue<TreeNode> preorderTraversal(TreeNode root) {
		Queue<TreeNode> ret = new LinkedList<TreeNode>();
		if(root == null) return(ret);
		Stack<TreeNode> s = new Stack<TreeNode>();
		s.push(root);
		while(s.size()>0) {
			TreeNode cur = s.pop();
			ret.offer(cur);
			if(cur.right!=null) s.push(cur.right);
			if(cur.left!=null) s.push(cur.left);
		}
		return(ret);
	}
	static boolean subtree2(TreeNode root, TreeNode subtree) {
		if (subtree == null) return(true);
		if (root == null) return(false);
		if(root.val == subtree.val) {
			if(matchTree(root, subtree)) return(true);
		}
		return(subtree2(root.left, subtree) || subtree2(root.right, subtree));
	}
	static boolean matchTree(TreeNode a, TreeNode b) {
		if(a == null && b != null) return(false);
		if(a != null && b == null) return(false);
		if(a == null && b == null) return(true);
		if(a.val == b.val) return(matchTree(a.left, b.left) && matchTree(a.right,b.right));
		else return(false);
	}

}
