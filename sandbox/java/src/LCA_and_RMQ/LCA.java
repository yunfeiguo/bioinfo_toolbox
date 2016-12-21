package LCA_and_RMQ;

import yunfeiImplementAlgs4.BinaryTreeNode;

/**
 * interface for lowest common ancestor problem.
 * @author guoy28
 *
 */
public interface LCA {
    /**
     * take two nodes in the tree, return their LCA,
     * i.e., parent farthest from root.
     * @param p
     * @param q
     * @return
     */
    public BinaryTreeNode query(BinaryTreeNode p, BinaryTreeNode q);
}
