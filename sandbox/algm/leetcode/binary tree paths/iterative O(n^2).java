/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode(int x) { val = x; }
 * }
 */
 import java.util.*;
public class Solution {
    public List<String> binaryTreePaths(TreeNode root) {
        List<String> ret = new LinkedList<String>();
        if (root == null) return(ret);
        Stack<TreeNode> s = new Stack<TreeNode>();
        Stack<String> curPath = new Stack<String>(); //record path from root to current node
        Map<TreeNode, Boolean> visited = new Hashtable<TreeNode, Boolean>();
        s.push(root);
        curPath.push(""+root.val);
        while(!s.isEmpty()) {
            TreeNode cur = s.peek();
            if(cur.left == null && cur.right == null) {
                //a leaf node
                ret.add(concatStr(curPath));
                curPath.pop();//pop out the leaf node from current paths
                s.pop();//pop out the leaf from nodes to be visited
            } else {
            //not a leaf node
            if (cur.left != null) {
                if (!(visited.containsKey(cur.left) && visited.get(cur.left))) {
                                 curPath.push(""+cur.left.val);
                s.push(cur.left);
                continue;   
                }
            }
            if (cur.right != null) {
                if (!(visited.containsKey(cur.right) && visited.get(cur.right))) {
                curPath.push(""+cur.right.val);
                s.push(cur.right);
                continue;
                }
            }
            curPath.pop();
            s.pop();
            }
            visited.put(cur,true);//this node has been visited
        }
        return(ret);
    }
    private String concatStr(Stack<String> s) {
        String ret = new String();
        int n = s.size();
        for (String i:s) {
        	if (n == 1) {
        		ret = ret + i;
        	} else {
        		ret = ret + i + "->";
        	}
        	n--;
        }
        return(ret);
    }
}
