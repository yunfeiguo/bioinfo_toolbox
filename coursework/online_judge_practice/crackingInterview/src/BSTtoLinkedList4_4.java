import java.util.LinkedList;
import java.util.Queue;
public class BSTtoLinkedList4_4 {
	public static LinkedList<LinkedList<Integer>> bstToLinkedList(TreeNode root) {
		LinkedList<LinkedList<Integer>> ret = new LinkedList<LinkedList<Integer>>(); 
		if(root == null) return(ret);
		Queue<TreeNode> q = new LinkedList<TreeNode>();
		q.offer(root);
		while(!q.isEmpty()) {
			int n = q.size();
			LinkedList<Integer> lvl = new LinkedList<Integer>();
			for(int i = 0; i < n; i++) {
				TreeNode cur = q.poll();
				lvl.add(cur.val);
				if(cur.left != null) q.offer(cur.left);
				if(cur.right != null) q.offer(cur.right);
			}
			ret.add(lvl);			
		}
		return(ret);
	}
	public static void main(String[] args) {
		TreeNode root = new TreeNode(3);
		root.left = new TreeNode(2);
		root.right = new TreeNode(4);
		LinkedList<LinkedList<Integer>> l = BSTtoLinkedList4_4.bstToLinkedList(root);
		for(LinkedList<Integer> i : l) {
			for (int j : i) {
				System.out.print(j + " ");
			}
			System.out.print("\n");
		}		
	}

}
