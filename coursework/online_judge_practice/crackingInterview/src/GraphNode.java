import java.util.List;
import java.util.LinkedList;
public class GraphNode {
	int val;
	List<GraphNode> neighbor = new LinkedList<GraphNode>();
	GraphNode (int x){
		this.val = x;
	}

}
