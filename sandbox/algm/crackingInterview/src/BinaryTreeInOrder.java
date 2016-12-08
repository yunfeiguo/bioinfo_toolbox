import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
public class BinaryTreeInOrder {
	/**
	 * recursive in-order traversal
	 * @param root
	 * @return
	 */
	public static List<Integer> inorderTraversal(TreeNode root) { 
		List<Integer> ret = new LinkedList<Integer>();
		if (root != null) {
			List<Integer> left = inorderTraversal(root.left);
			List<Integer> right = inorderTraversal(root.right);
			if (left != null) ret.addAll(left);
			ret.add(root.val);
			if (right != null) ret.addAll(right);       
		}
		return(ret);
	}
	/**
	 * iterative in-order traversal
	 * @param root
	 * @return
	 */
	public static List<Integer> inorderTraversalIter(TreeNode root) {
		List<Integer> ret = new LinkedList<Integer>();		
		Stack<TreeNode> s = new Stack<TreeNode>();
		TreeNode cur = root;
		while (cur != null) {
			s.push(cur);
			cur = cur.left;
			if(cur != null) {
				continue;
			} else {
				while (!s.isEmpty()) {
					TreeNode toBeShown = s.pop();
					ret.add(toBeShown.val);
					if (toBeShown.right != null) {
						cur = toBeShown.right;
						break;
					}
				}
			}
		}
		return(ret);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(0);
		root.right = new TreeNode(2);
		root.right.left = new TreeNode(3);
		List<Integer> traversal = BinaryTreeInOrder.inorderTraversalIter(root);
		for (int i:traversal) {
			System.out.println(i);
		}
	}
}
