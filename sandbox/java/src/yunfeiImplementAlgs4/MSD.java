package yunfeiImplementAlgs4;

public class MSD {
	private static int R = 65535; //number of unique chars
	private static final int M = 15; //cutoff for direct sorting
	private static String[] aux;
	
	private static int charAt(String s, int d) {
		if (d < s.length()) {
			return(s.charAt(d));
		} else {
			return(-1);
		}
	}
	public static void sort(String[] a) {
		int N = a.length;
		aux = new String[N];
		sort(a,0,N-1,0);
	}
	private static void sort(String[] a, int lo, int hi, int d) {
		int[] count = new int[R+2];
		if (lo + M >= hi) {
			insertion(a, lo, hi, d);
			return;
		}
		//we may use interstion sorting to speed up	
		//get frequency
		for(int i = lo; i <= hi; i++) {
			count[charAt(a[i], d)+2]++;
		}
		//get starting position
		for (int i = 0; i < R + 1; i++) {
			count[i+1] += count[i];
		}
		//distribute elements, convert counts to index
		for (int i = lo; i <= hi; i++) {
			aux[count[charAt(a[i],d) + 1]++] = a[i];
		}
		//copy back
		for (int i = lo; i <= hi; i++)	{
			a[i] = aux[i - lo];
		}
		
		//do this for every character
		for (int i = 0; i < R; i++) {
			sort(a,count[i] + lo,count[i+1]-1 + lo,d+1);
		}	
	}
    // insertion sort a[lo..hi], starting at dth character
	private static void insertion(String[] a, int lo, int hi, int d	) {
		for (int i = lo; i <= hi; i++) {
			for (int j = i; j > lo && less(a[j],a[j-1],d); j--) {
				exch(a,j,j-1);
			}
		}
	}
	/* exchange elements in a string array
	 * @param a string array
	 * @param i first index
	 * @param j second index
	 */
	private static void exch(String[] a, int i, int j) {
		String tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	}
	private static Boolean less(String v, String w, int idx) {
		if (idx < v.length() && idx < w.length()) {
			return(v.substring(idx).compareTo(w.substring(idx)) < 0);
		} else if (idx < v.length()) {
			return(false);
		} else if (idx < w.length()) {
			return(true);
		} else {
			return(false);
		}
	}
}
