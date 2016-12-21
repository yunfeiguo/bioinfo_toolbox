package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.Quick3string;
import edu.princeton.cs.algs4.MSD;
import edu.princeton.cs.algs4.StdIn;

public class CompareAlgs {
	public static void main(String[] args) {		
		String[] a = StdIn.readAllStrings();
		
		long start = System.currentTimeMillis();
		Quick3string.sort(a);
		long end = System.currentTimeMillis();
		System.out.println(end - start);
		
		/*
		long start = System.currentTimeMillis();
		MSD.sort(a);
		long end = System.currentTimeMillis();
		System.out.println(end - start);
		*/
	}
}
