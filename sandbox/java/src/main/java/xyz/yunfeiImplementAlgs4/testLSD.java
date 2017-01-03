package yunfeiImplementAlgs4;
import edu.princeton.cs.algs4.*;
public class testLSD {
	public static void main(String[] args) {
		String[] a = new String[5];
		a[0] = "4PGC938";
		a[1] = "2IYE230";
		a[2] = "3CI0720";
		a[3] = "3CI0710";
		a[4] = "10HV845";
		LSD.sort(a, 7);
		for (int i = 0; i < a.length; i++) {
			StdOut.println(a[i]);
		}
	}

}
