import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
public class BinaryTreePreOrder {
	/**
	 * recursive preorder binary tree traversal
	 * @param root
	 * @return
	 */
	public static List<Integer> preorderTraversal(TreeNode root) {
		List<Integer> ret = new LinkedList<Integer>();
		if(root != null) {
			ret.add(root.val);
			ret.addAll(preorderTraversal(root.left));
			ret.addAll(preorderTraversal(root.right));
		}
		return(ret);
	}
	/**
	 * iterative pre-order traversal
	 * @param root
	 * @return
	 */
	public static List<Integer> preorderTraversalIter(TreeNode root) {
		List<Integer> ret = new LinkedList<Integer>();
		if(root == null) return(ret);
		Stack<TreeNode> s = new Stack<TreeNode>();
		s.push(root);
		while (!s.isEmpty()) {
			TreeNode cur = s.pop();
			if (cur.right != null) s.push(cur.right);
			if (cur.left != null) s.push(cur.left);
			ret.add(cur.val);
		}			
		return(ret);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(1);
		root.right = new TreeNode(2);
		root.right.left = new TreeNode(3);
		List<Integer> l = BinaryTreePreOrder.preorderTraversalIter(root);
		for(int i:l) {
			System.out.println(i);
		}
	}
}
