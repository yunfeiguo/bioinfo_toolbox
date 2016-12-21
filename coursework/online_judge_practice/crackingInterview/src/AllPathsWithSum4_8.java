import java.util.*;
public class AllPathsWithSum4_8 {
	public static void main(String[] args) {
		TreeNode root = new TreeNode(5);
		root.left = new TreeNode(-1);	
		root.left.right = new TreeNode(1);	
		root.left.left = new TreeNode(3);
		root.left.left.left = new TreeNode(-1);
		root.left.left.left.left = new TreeNode(1);
		root.right = new TreeNode(1);		
		root.right.right = new TreeNode(7);
		LinkedList<TreeNode> prevNodes = new LinkedList<TreeNode>();
		AllPathsWithSum4_8.findSum(root,7,prevNodes, 0);
		/* ans
		 * 5 -1 3 
		 * 5 -1 3 -1 1 
		 * 7 
		 */
	}
	/**
	 * print all paths with sum t starting from level l
	 * @param root
	 * @param t
	 * @param prev: store previous nodes on the path
	 * @param l
	 */
	public static void findSum(TreeNode root, int t, LinkedList<TreeNode> prev, int l) {
		if (root == null) return;
		int sum = 0;
		prev.add(root);
		for(int i = l; i > -1; i--) {
			sum = sum + prev.get(i).val;
			if(sum == t) {				
				print(prev, i, l);
			}
		}	
		findSum(root.left, t, (LinkedList<TreeNode>) prev.clone(), l + 1);
		findSum(root.right, t, (LinkedList<TreeNode>) prev.clone(), l + 1);
		return;
	}
	static void print(LinkedList<TreeNode> l, int s, int e) {
		for(int i = s; i <= e; i++) {
			System.out.print(l.get(i).val + " ");
		}
		System.out.println();
	}
}
