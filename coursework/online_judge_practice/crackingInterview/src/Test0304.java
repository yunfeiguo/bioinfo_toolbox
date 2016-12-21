
public class Test0304 {
    public int heightIterative() {
        Queue q = new LinkedList();
        BTreeNode node;
        int currentLevelNodes, nextLevelNodes,depth;
        if (root == null)
            return 0;
        q.add(root);
        currentLevelNodes = 1;depth=0;
        nextLevelNodes = 0;
        while (!q.isEmpty()) {
            node = q.remove();
            currentLevelNodes–;
            if (node.left != null) {
                q.add(node.left);
                nextLevelNodes++;
            }
            if (node.right != null) {
                q.add(node.right);
                nextLevelNodes++;
            }
            if (currentLevelNodes == 0) {
                depth++;
                currentLevelNodes = nextLevelNodes;
                nextLevelNodes = 0;
            }
        }
        return depth;
    }
}
