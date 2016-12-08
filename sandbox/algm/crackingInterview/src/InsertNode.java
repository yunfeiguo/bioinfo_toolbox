import java.util.List;
public class InsertNode {
	/**
	 * recursive solution
	 * assume no equality
	 */
	public static TreeNode insertNode(TreeNode root, TreeNode t) {
		if (root == null) {
			return(t);
		}
		if (t == null) {
			return(root);
		}
		if (root.val > t.val) {
			root.left = insertNode(root.left,t);			
		} else {
			root.right = insertNode(root.right,t);
		}
		return(root);		
	}
	public static TreeNode insertNodeIter(TreeNode root, TreeNode t) {
		if (root == null) {
			return(t);
		}
		if (t == null) {
			return(root);
		}
		TreeNode prev = root;
		TreeNode cur = root;
		while (cur != null) {
			prev = cur;
			if (t.val > cur.val) {				
				cur = cur.right;
			} else {				
				cur = cur.left;
			}			
		}
		if (t.val > prev.val) {
			prev.right = t;
		} else {
			prev.left = t;
		}
		return(root);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(10);
		root.right = new TreeNode(20);
		InsertNode.insertNodeIter(root,new TreeNode(15));
		List<Integer> l = BinaryTreeInOrder.inorderTraversalIter(root);
		for(int i:l){
			System.out.println(i);
		}
	}
	
}
