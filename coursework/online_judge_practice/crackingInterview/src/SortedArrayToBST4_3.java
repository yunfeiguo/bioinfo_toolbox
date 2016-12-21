import java.util.LinkedList;
import java.util.Stack;

public class SortedArrayToBST4_3 {
	/**
	 * take an ascending integer array, return root node of a BST
	 * @param a
	 * @return
	 */
	public static TreeNode bst (int[] a, int lo, int hi) {
		if(lo > hi) {
			return(null);
		}
		if (lo == hi) {
			return(new TreeNode(a[lo]));
		}
		int mid = lo + (hi - lo)/2;
		TreeNode root = new TreeNode(a[mid]);
		root.left = bst(a, lo, mid - 1);
		root.right = bst(a, mid + 1, hi	);
		return(root);
	}
	/**
	 * use iterative method for bst construction
	 * @param a
	 * @return
	 */
	public static TreeNode bstIter(int[] a) {
		if(a.length == 0) return(null);
		TreeNode root = new TreeNode(0);
		Stack<TreeNode> nodeStack = new Stack<TreeNode>();
		Stack<Integer> leftIdx = new Stack<Integer>();
		Stack<Integer> rightIdx = new Stack<Integer>();
		nodeStack.push(root);
		leftIdx.push(0);
		rightIdx.push(a.length - 1);
		while(!nodeStack.isEmpty()) {
			TreeNode cur = nodeStack.pop();
			int l = leftIdx.pop();
			int r = rightIdx.pop();
			int mid = l + (r -l)/2;
			cur.val = a[mid];
			
			if(l <= mid - 1) {
				cur.left = new TreeNode(0);				 
				nodeStack.push(cur.left);
				leftIdx.push(l);
				rightIdx.push(mid - 1);
			}
			if(r >= mid + 1) {
				cur.right = new TreeNode(0);
				nodeStack.push(cur.right);
				leftIdx.push(mid + 1);
				rightIdx.push(r);
			}
		}
		return(root);
	}
	public static void main(String[] args) {
		int[] a = {1,2,3,4,5};
		//TreeNode t = SortedArrayToBST.bst(a, 0, a.length-1);
		TreeNode t = SortedArrayToBST4_3.bstIter(a);
		if (t != null) {
			t.printTree();
		}
	}
}
