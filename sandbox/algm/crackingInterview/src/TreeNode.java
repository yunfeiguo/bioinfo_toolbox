import java.util.Queue;
import java.util.LinkedList;
public class TreeNode {
	int val;
	TreeNode left;
	TreeNode right;
	TreeNode parent;
	TreeNode (int x) { val = x;}
	public void printTree () {
		//pre-order traversal
		Queue<TreeNode> q = new LinkedList<TreeNode>();		
		q.offer(this);
		TreeNode cur;
		int n;
		while(!q.isEmpty()) {
			n = q.size();
			for (int i = 0; i < n; i++) {
				cur = q.poll();
				if (cur == null) {
					System.out.print("null");
				} else {				
					q.offer(cur.left);
					q.offer(cur.right);				
					System.out.print(cur.val);
				}
				System.out.print(" ");
			}	
		}
		System.out.print("\n");
		return;
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(0);
		root.printTree();
	}
}
