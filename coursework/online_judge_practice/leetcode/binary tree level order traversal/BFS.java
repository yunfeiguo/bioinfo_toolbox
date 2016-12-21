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
}
