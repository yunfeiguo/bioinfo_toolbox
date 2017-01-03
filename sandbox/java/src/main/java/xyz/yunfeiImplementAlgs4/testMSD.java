package yunfeiImplementAlgs4;
import edu.princeton.cs.algs4.*;
public class testMSD {
	public static void main(String[] args) {
		String[] a = new String[5];
		a[0] = "4PGC93";
		a[1] = "2IYE2";
		a[2] = "3CI0";
		a[3] = "3CI0710";
		a[4] = "10HV845";
		MSD.sort(a);
		for (int i = 0; i < a.length; i++) {
			StdOut.println(a[i]);
		}
	}
}
