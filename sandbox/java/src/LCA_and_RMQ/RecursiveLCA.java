package LCA_and_RMQ;

import yunfeiImplementAlgs4.BinaryTreeNode;

/*
 * recursive method for LCA problem
 */
public class RecursiveLCA implements LCA {
    private BinaryTreeNode root;
    public RecursiveLCA(BinaryTreeNode root) {
        this.root = root;
    }
    /**
     * given root and p, q in the tree, find lowest common ancestor of p and q
     * @param root
     * @param p
     * @param q
     * @return
     */
    public BinaryTreeNode query(BinaryTreeNode p, BinaryTreeNode q) {
        return query(root, p, q);
    }
    private BinaryTreeNode query(BinaryTreeNode root, BinaryTreeNode p, BinaryTreeNode q) {       
        if (root == null || root == p || root == q) {
            return root;
        }
        BinaryTreeNode leftAncestor = query(root.left, p, q);
        BinaryTreeNode rightAncestor = query(root.right, p, q);
        if (leftAncestor == null) {
            return rightAncestor;
        } else if (rightAncestor == null){
            return leftAncestor;
        } else {
            return root;
        }
    }
}
