
public class FindNextNodeInBST4_5 {
	public static TreeNode getNext(TreeNode n) {
		if(n == null) return(n);
		if(n.parent == null) {
			//this is the root node, return leftmost node in right subtree
			n = n.right;
			while(n != null && n.left != null) {
				n = n.left;
			}
			return(n);
		}
		//if this is not a root node
		//search in its parent, until no parent has been found (encounter the root)
		//or parent.val > cur.val
		if(n.parent != null) {
			int t = n.val;
			while (n.parent != null && n.val <= t) {
				n = n.parent;
			}
			if (n.parent == null && n.val <= t) {
				n = null; //we get to the root, but still no valid successor
			}
		}
		return(n);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(5);
		root.left = new TreeNode(3);
		root.left.parent = root;
		root.left.right = new TreeNode(4);
		root.left.right.parent = root.left;
		root.left.left = new TreeNode(1);
		root.left.left.parent = root.left;
		root.right = new TreeNode(7);
		root.right.parent = root;
		root.right.left = new TreeNode(6);
		root.right.left.parent = root.right;
		root.right.right = new TreeNode(9);
		root.right.right.parent = root.right;
		assert(FindNextNodeInBST4_5.getNext(root).val == 6);	
		assert(FindNextNodeInBST4_5.getNext(root.left.right).val == 5);
		assert(FindNextNodeInBST4_5.getNext(root.right.right) == null);
	}

}
