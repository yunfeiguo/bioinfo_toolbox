package LCA_and_RMQ;

import edu.princeton.cs.algs4.In;
import yunfeiImplementAlgs4.BinaryTreeNode;

public class TestLCA {
    public static void main(String[] args) {
        for (int i = 0; i < args.length; i++) {   
            In in = new In(args[i]);
            BinaryTreeNode root = BinaryTreeNode.deserialize(in.readAll().trim());
            BinaryTreeNode leftLeaf = root;
            while(leftLeaf.left != null) {
                leftLeaf = leftLeaf.left;
            }
            BinaryTreeNode rightLeaf = root;
            while(rightLeaf.right != null) {
                rightLeaf = rightLeaf.right;
            }
            //LinearLCA llca = new LinearLCA(root);
            RecursiveLCA rlca = new RecursiveLCA(root);
            //System.out.println(llca.query(root, leftLeaf) == root);
            System.out.println(rlca.query(root, leftLeaf) == root);
            //System.out.println(llca.query(root, rightLeaf) == root);
            System.out.println(rlca.query(root, rightLeaf) == root);
            //System.out.println(llca.query(rightLeaf, leftLeaf) == root);
            System.out.println(rlca.query(rightLeaf, leftLeaf) == root);
        }
    }
}
