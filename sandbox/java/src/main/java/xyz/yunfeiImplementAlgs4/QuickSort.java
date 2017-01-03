package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdRandom;

public class QuickSort {
	public static void sort(Comparable[] a) {

		StdRandom.shuffle(a);
		sort(a,0,a.length-1);
	}
	private static void printA(Comparable[] a, int lo, int hi) {
		System.out.println("sorting:");
		for (int i = lo; i <= hi; i++) {
			System.out.println(a[i]);
		}
	}
	private static void sort(Comparable[] a, int lo, int hi) {
		printA(a, lo, hi);
		if (lo >= hi) {
			return;
		}
		int j = partition(a,lo,hi);
		sort(a, lo, j - 1);
		sort(a, j + 1, hi);
	}
	private static int partition(Comparable[] a, int lo, int hi) {
		Comparable pivot = a[hi];
		int i = lo;
		for (int j = lo; j < hi; j++) {
			if (less(a[j],pivot)) {
				exch(a,i++,j);				
			}
		}
		exch(a, hi, i);
		return(i);		
	}
	private static Boolean less(Comparable a, Comparable b) {
		return(a.compareTo(b) < 0);
	}
	private static void exch(Comparable[] a, int i, int j) {
		Comparable tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	}
}
