import java.util.*;
public class FindNearestCommonAncestor4_6 {
	public static void main(String[] args) {
		TreeNode root = new TreeNode(5);
		root.left = new TreeNode(3);	
		root.left.right = new TreeNode(4);	
		root.left.left = new TreeNode(1);	
		root.right = new TreeNode(7);	
		root.right.left = new TreeNode(6);	
		root.right.right = new TreeNode(9);
		System.out.println(FindNearestCommonAncestor4_6.findComAncestor(root,root,root.right).val == root.val);			
		System.out.println(FindNearestCommonAncestor4_6.findComAncestor(root,root.left.right,root.right.left).val == root.val);
		System.out.println(FindNearestCommonAncestor4_6.findComAncestor(root,root.left,root.left.right).val == root.left.val);
		System.out.println(FindNearestCommonAncestor4_6.findComAncestor(root,root.left.left,root.right.right).val == root.val);
	}
	public static TreeNode findComAncestor(TreeNode root, TreeNode a, TreeNode b) {		
		if (root == null) return(null);
		if (a == b) return(a);
		if (a == null) return(b);
		if (b == null) return(a);
		if (root == a || root == b) return(root);
		boolean leftRootA = isChild(root.left, a);
		boolean leftRootB = isChild(root.left, b);
		if (leftRootA == leftRootB) {
			if (leftRootA == false) {
				return(findComAncestor(root.right, a, b));
			} else {
				return(findComAncestor(root.left, a, b));	
			}
		} else {
			return(root);
		}
	}
	/**
	 * determine if x is a child(grandchild) of root
	 * @param root
	 * @param x
	 * @return
	 */
	private static boolean isChild(TreeNode root, TreeNode x) {
		//use BFS to traverse
		Queue<TreeNode> q = new LinkedList<TreeNode>();
		if (root == null || x == null) return(false);
		q.offer(root);
		while(!q.isEmpty()) {
			int n = q.size();
			for(int i = 0; i < n; i++) {
				TreeNode cur = q.poll();
				if (cur == x) return(true);
				if(cur.left != null) q.offer(cur.left);
				if(cur.right != null) q.offer(cur.right);
			}
		}
		return(false);		
	}
}
