package yunfeiImplementAlgs4;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;

import edu.princeton.cs.algs4.In;

public class BinaryTreeNode {
    public BinaryTreeNode left;
    public BinaryTreeNode right;    
    public int key;
    public BinaryTreeNode(Integer key) {
        this.key = key;
        this.left = null;
        this.right = null;
    }        
    /**
     * return Euler Tour in an array. An Euler Tour
     * traverses every edge of a directed graph
     * exactly once and records which nodes are visited.
     * each undirected edge consists of two directed 
     * edges of opposite directions.
     * 
     * @param root
     * @return
     */
    public static BinaryTreeNode[] EulerTour(BinaryTreeNode root) {
        ArrayList<BinaryTreeNode> tour = new ArrayList<BinaryTreeNode>();
        EulerTour(tour, root);
        return tour.toArray(new BinaryTreeNode[0]);
    }
    private static void EulerTour(ArrayList<BinaryTreeNode> tour, BinaryTreeNode root) {
        if (root == null) {
            return;
        }
        tour.add(root);
        EulerTour(tour, root.left);
        if (root.left != null) {
            tour.add(root);
        }
        EulerTour(tour, root.right);
        if (root.right != null) {
            tour.add(root);
        }
    }
    /**
     * output levels (distance from root) in Euler tour order
     * @return
     */
    public static int[] getLevels(BinaryTreeNode root, BinaryTreeNode[] euler) {
        if (root == null || euler == null || euler.length == 0) {
            throw new IllegalArgumentException();
        }
        Map<BinaryTreeNode, Integer> node2levels = BinaryTreeNode.getLevels(root);
        int[] levels = new int[euler.length];
        for (int i = 0; i < euler.length; i++) {
            levels[i] = node2levels.get(euler[i]);            
        }
        return levels;
    }
    /**
     * output levels in a Map with key as nodes
     * and values as levels. Use BFS (level-order
     * traversal) to find all levels.
     */
    public static Map<BinaryTreeNode, Integer> getLevels(BinaryTreeNode root) {
        Queue<BinaryTreeNode> toVisit = new LinkedList<BinaryTreeNode>();
        Map<BinaryTreeNode, Integer> node2level = new HashMap<BinaryTreeNode, Integer>();
        int currentLevel = 0;
        toVisit.offer(root);
        while(!toVisit.isEmpty()) {
            int currentLevelSize = toVisit.size();
            for (int i = 0; i < currentLevelSize; i++) {
                BinaryTreeNode current = toVisit.poll();
                if(current.left != null) {
                    toVisit.add(current.left);
                }
                if (current.right != null) {
                    toVisit.add(current.right);
                }
                node2level.put(current, currentLevel);
            }
            currentLevel++;
        }
        return node2level;
    }
    /**
     * take root of a binary tree and serialize it into text
     * @param root
     * @return
     */
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
    /**
     * simple unit tests
     * @param args
     */
    public static void main(String[] args) {
        /*test serialization deserialization*/
        for (int i = 0; i < args.length; i++) {
            System.out.println("processing " + args[i]);
            In in = new In(args[i]);
            System.out.println(serialize(deserialize(in.readAll().trim())));   
        }        
        System.out.println(serialize(deserialize("1,2,3")));
        System.out.println(serialize(deserialize("1,2,#,3,#")));
        /*test Euler tour*/
        BinaryTreeNode test = new BinaryTreeNode(0);
        test.left = new BinaryTreeNode(1);
        test.right = new BinaryTreeNode(2);
        BinaryTreeNode[] eulerTour = BinaryTreeNode.EulerTour(test);
        System.out.println(eulerTour[0] == test);
        System.out.println(eulerTour[1] == test.left);
        System.out.println(eulerTour[2] == test);
        System.out.println(eulerTour[3] == test.right);
        System.out.println(eulerTour[4] == test);
        /*test levels*/
        Map<BinaryTreeNode, Integer> node2level = BinaryTreeNode.getLevels(test);
        System.out.println(node2level.get(test) == 0);
        System.out.println(node2level.get(test.left) == 1);
        System.out.println(node2level.get(test.right) == 1);
    }
}
