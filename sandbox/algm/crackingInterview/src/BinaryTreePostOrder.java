import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

public class BinaryTreePostOrder {
	/**
	 * recursive postorder binary tree traversal
	 * @param root
	 * @return
	 */
	public static List<Integer> postorderTraversal(TreeNode root) {
		List<Integer> ret = new LinkedList<Integer>();
		if(root != null) {
			ret.addAll(postorderTraversal(root.left));
			ret.addAll(postorderTraversal(root.right));
			ret.add(root.val);
		}
		return(ret);
	}
	/**
	 * iterative post-order traversal
	 * we just get root-right-left and then reverse it!
	 * @param root
	 * @return
	 */
	public static List<Integer> postorderTraversalIter(TreeNode root) {
		List<Integer> ret = new LinkedList<Integer>();
		if(root == null) return(ret);
		Stack<TreeNode> s = new Stack<TreeNode>();
		s.push(root);
		while (!s.isEmpty()) {
			TreeNode cur = s.pop();
			if (cur.left != null) s.push(cur.left);
			if (cur.right != null) s.push(cur.right);
			ret.add(cur.val);
		}		
		Collections.reverse(ret);
		return(ret);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(1);
		root.right = new TreeNode(2);
		root.right.left = new TreeNode(3);
		List<Integer> l = BinaryTreePostOrder.postorderTraversalIter(root);
		for(int i:l) {
			System.out.println(i);
		}
	}
}
