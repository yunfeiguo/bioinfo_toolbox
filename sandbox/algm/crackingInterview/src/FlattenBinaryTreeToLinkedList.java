import java.util.Stack;
public class FlattenBinaryTreeToLinkedList {
	/**
	 * use pre-order traversal to flatten binary tree
	 * @param root
	 */
	public static void flatten ( TreeNode root) {
		if(root == null) return;
		Stack<TreeNode> s = new Stack<TreeNode>();
		TreeNode prev = root;
		s.push(root);
		while(!s.isEmpty()) {
			TreeNode cur = s.pop();
			if(cur.right != null) s.push(cur.right);
			if(cur.left != null) s.push(cur.left);
			if(cur != prev) prev.right = cur;
			prev = cur;
			prev.left = null;
			prev.right = null;
		}
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(1);
		root.left = new TreeNode(2);
		root.left.left = new TreeNode(3);
		root.left.right = new TreeNode(4);
		root.right = new TreeNode(5);
		root.right.right = new TreeNode(6);
		FlattenBinaryTreeToLinkedList.flatten(root);
		root.printTree();
	}

}
