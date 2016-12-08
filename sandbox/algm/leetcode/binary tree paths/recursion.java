/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode(int x) { val = x; }
 * }
 */
//time complexity: nlogn
//space: O(n) without stack space
 import java.util.*;
public class Solution {
    public List<String> binaryTreePaths(TreeNode root) {
        List<String> ret = new LinkedList<String>();
        if (root == null) return(ret);
        if (root.left == null && root.right == null) {
            ret.add(""+root.val);
        }
        
        for (String i:binaryTreePaths(root.left)) {
            ret.add(root.val+"->"+i);
        
        }
        
        for (String i:binaryTreePaths(root.right)) {
            ret.add(root.val+"->"+i);
        
        }
        return(ret);
    }
}
