package yunfeiImplementAlgs4;
import java.util.*;
public class Tree2Text {   
    public static String serialize(BinaryTreeNode root) {
        StringBuilder sb = new StringBuilder();
        java.util.Queue<BinaryTreeNode> queue = new LinkedList<BinaryTreeNode>();
        queue.offer(root);
        
        while(!queue.isEmpty()){
            BinaryTreeNode top = queue.poll();
            if(top == null){
                sb.append("#,");
            }
            else{
                sb.append(top.key + ",");
                queue.offer(top.left);
                queue.offer(top.right);
            }
        }
        
        return removeTrailingNulls(sb.toString());
    }
    /**
     * remove trailing # signs and commas 
     */
    private static String removeTrailingNulls(String s) {
        int newLength = s.length();
        for (int i = s.length() - 1; i >= 0; i--) {
            if (s.charAt(i) == '#' || s.charAt(i) == ',') {
                newLength = i;
            } else {
                break;
            }
        }
        return s.substring(0, newLength);
    }
    // Decodes your encoded data to tree.
    public static BinaryTreeNode deserialize(String data) {
        if(data.isEmpty() || data.charAt(0) == '#'){
            return null;
        }

        String[] arr = data.split(",");
        int idx = 0;

        BinaryTreeNode root = new BinaryTreeNode(Integer.valueOf(arr[0]));
        java.util.Queue<BinaryTreeNode> queue = new LinkedList<BinaryTreeNode>();
        queue.offer(root);

        while(!queue.isEmpty()){
            BinaryTreeNode top = queue.poll();
            if(top == null){
                continue;
            }

            if (idx < arr.length - 1 && !arr[++idx].equals("#")) {
                top.left = new BinaryTreeNode(Integer.valueOf(arr[idx]));
            }
            if (idx < arr.length - 1 && !arr[++idx].equals("#")) {
                top.right = new BinaryTreeNode(Integer.valueOf(arr[idx]));
            }            
            queue.offer(top.left);
            queue.offer(top.right);
        }

        return root;
    }
    public static void main(String[] args) {
        System.out.println(Tree2Text.serialize(Tree2Text.deserialize("1,2,3")));
        System.out.println(Tree2Text.serialize(Tree2Text.deserialize("1,2,#,3,#")));
    }
}