package yunfeiImplementAlgs4;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import java.util.List;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Queue;

public class TestInstanceArray { 
    //private static final int x = nonStaticNum();
    private static int staticNum(){
        return 1;
    }
    private int nonStaticNum() {
        return staticNum();
    }
    private static int staticUsingNonStatic() {
        //return nonStaticNum();
        return staticNum();
    }
    public static void main(String[] args) {
        Deque<Integer> x = new LinkedList<Integer> ();
        x.offerLast(1); //same as offer() in Queue
        x.offerLast(2);
        System.out.println(1 == x.pollFirst()); //same as poll() in Queue
        System.out.println(2 == x.peekFirst()); //same as peek() in Queue

        x.addFirst(-1); //push() in Stack
        System.out.println(x.removeFirst() == -1); //same as pop() in Stack
        




        
    }

}
