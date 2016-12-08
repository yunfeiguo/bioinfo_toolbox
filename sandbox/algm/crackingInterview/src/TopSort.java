import java.util.ArrayList;
public class TopSort {
    public ArrayList<DirectedGraphNode> topSort(ArrayList<DirectedGraphNode> graph) {
        ArrayList<DirectedGraphNode> ret = new ArrayList<DirectedGraphNode> ();
        if (graph == null || graph.size() == 0) {
            return ret;
        }        
        //assume no cycle exists
        topSort(ret, graph);
        return ret;
    }
    private void topSort(ArrayList<DirectedGraphNode> ret, ArrayList<DirectedGraphNode> g) {
        if (g.size() == 1) {
            //add the root node
            ret.add(0, g.get(0));
            return;
        }
        //add all nodes with no children to head of ret
        //then remove these nodes from g, and also remove them from the neighbors
        ArrayList<DirectedGraphNode> newG = new ArrayList<DirectedGraphNode> ();
        for (DirectedGraphNode n : g) {
            if (n.neighbors == null || n.neighbors.size() == 0) {
                ret.add(0, n);
            } else {
                removeHangs(n);
                newG.add(n);
            }
        }
        g = newG;
        topSort(ret, g);
    }
    /*not work*/
    private void removeHangs(DirectedGraphNode n) {
        ArrayList<Integer> toRemove = new ArrayList<Integer> ();
        //remove nodes with zero neighbors in neighbors
        for (int i = n.neighbors.size() - 1; i >= 0; i--) {
            DirectedGraphNode cur = n.neighbors.get(i);
            if (cur.neighbors == null || cur.neighbors.size() == 0) {
                toRemove.add(i);
            }
        }
        for (int i : toRemove) {
            
        }
    }
    public static void main(String[] args) {
        DirectedGraphNode one = new DirectedGraphNode(1);
        DirectedGraphNode two = new DirectedGraphNode(2);
        DirectedGraphNode three = new DirectedGraphNode(3);
        DirectedGraphNode four = new DirectedGraphNode(4);
        one.neighbors.add(two);
        one.neighbors.add(three);
        two.neighbors.add(three);
        two.neighbors.add(four);
        three.neighbors.add(four);
        ArrayList<DirectedGraphNode> all = new ArrayList<DirectedGraphNode> ();
        all.add(four);
        all.add(two);
        all.add(three);
        all.add(one);
        TopSort test = new TopSort();
        for (DirectedGraphNode g : test.topSort(all)) {
            System.out.println(g.label);
        }               
    }
}
