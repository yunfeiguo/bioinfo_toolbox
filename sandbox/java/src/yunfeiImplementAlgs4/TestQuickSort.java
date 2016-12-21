package yunfeiImplementAlgs4;
import edu.princeton.cs.algs4.*;
public class TestQuickSort {
	public static void main(String[] args) {
		Integer[] a = new Integer[5];
		a[0] = 1;
		a[1] = 1;
		a[2] = 2;
		a[3] = 1;
		a[4] = 1;
		QuickSort.sort(a);	
		System.out.println("result:");
		for (int i = 0; i < a.length; i++) {
			StdOut.println(a[i]);
		}
	}

}
