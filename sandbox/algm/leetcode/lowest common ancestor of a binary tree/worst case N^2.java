/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode(int x) { val = x; }
 * }
 */
public class Solution {
    public TreeNode lowestCommonAncestor(TreeNode root, TreeNode a, TreeNode b) {
        		if (root == null) return(null);
		if (a == b) return(a);
		if (a == null) return(b);
		if (b == null) return(a);
		if (root == a || root == b) return(root);
		boolean leftRootA = isChild(root.left, a);
		boolean leftRootB = isChild(root.left, b);
		if (leftRootA == leftRootB) {
			if (leftRootA == false) {
				return(lowestCommonAncestor(root.right, a, b));
			} else {
				return(lowestCommonAncestor(root.left, a, b));	
			}
		} else {
			return(root);
		}
    }
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
