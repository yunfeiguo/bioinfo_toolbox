package LCA_and_RMQ;

import yunfeiImplementAlgs4.BinaryTreeNode;
import java.util.Map;
import java.util.HashMap;
/**
 * implementation of Michael Bender and Martin Farach-Colton's algorithm
 *
 * given a binary tree, use O(n) time for precomputing, such that
 * any query for lowest common ancestor of two nodes u and v only 
 * costs O(1) time.
 * 
 * @author guoy28
 */
public class LinearLCA implements LCA {
    /* 
     * algorithm
     * 
     * first construct an array of Euler tour E[0..2n-1]    
     * then construct an array of representative nodes R[0..n-1]
     * (first index of a node in Euler tour). The array can be substituted
     * by map if tree nodes are not indexed by integers. 
     * then construct an array of levels (distance from root) L[0..2n-1]
     * for all elements in Euler tour. 
     * for each query (u,v), lowest common ancestor the node whose level is 
     * smallest in L[R[u],R[v]]. Therefore we convert LCA problem to
     * RMQ problem, which can be solved in <O(nlgn), O(1)> or even <O(n), 
     * O(1)> time.
     */
    private BinaryTreeNode[] euler;
    private Map<BinaryTreeNode, Integer> node2representatives;
    private int[] levels;
    private RMQ rmq;
    
    /**
     * initialize Euler tour, representatives, levels
     * and RMQ object
     * @param root
     */
    public LinearLCA(BinaryTreeNode root) {
        euler = BinaryTreeNode.EulerTour(root);
        node2representatives = getRepresentatives(euler);      
        levels = BinaryTreeNode.getLevels(root, euler);
        this.rmq = new NormalizedRMQ(this.levels);
    }
    private Map<BinaryTreeNode, Integer> getRepresentatives(BinaryTreeNode[] euler) {
        Map<BinaryTreeNode, Integer> representatives = new HashMap<BinaryTreeNode, Integer>();
        for (int i = 0; i < euler.length; i++) {
            if (!representatives.containsKey(euler[i])) {
                representatives.put(euler[i], i);
            }
        }
        return representatives;
    }
    public BinaryTreeNode query(BinaryTreeNode u, BinaryTreeNode v) {
        if (!node2representatives.containsKey(u) || !node2representatives.containsKey(v)) {
            throw new IllegalArgumentException();
        }
        int firstU = node2representatives.get(u);
        int firstV = node2representatives.get(v);
        return euler[rmq.min(Math.min(firstU, firstV), Math.max(firstU,firstV))];        
    }
    /**
     * simple unit tests
     * @param args
     */
    public static void main(String[] args) {
        BinaryTreeNode root = new BinaryTreeNode(0);
        root.left = new BinaryTreeNode(1);
        root.right = new BinaryTreeNode(2);
        root.right.right = new BinaryTreeNode(5);
        root.right.left = new BinaryTreeNode(3);
        root.right.left.right = new BinaryTreeNode(4);
        LinearLCA lca = new LinearLCA(root);
        System.out.println(lca.query(root.left, root.right) == root);
        System.out.println(lca.query(root.left, root) == root);
        System.out.println(lca.query(root.right.left.right, root.right.right) == root.right);
    }
}