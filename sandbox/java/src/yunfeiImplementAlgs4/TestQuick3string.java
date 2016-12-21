package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.Quick3string;

public class TestQuick3string {
	public static void main(String[] args) {
		String[] a = {"ab","abd","abc"};
		//Quick3String2.sort(a);
		Quick3string.sort(a);
		System.out.println(String.join(" ", a));
	}
}
