package yunfeiImplementAlgs4;

public class Quick3String2 {
	/**
	 * sort strings using 3-way quicksort and MSD
	 * @param s
	 */
	public static void sort(String[] s) {
		sort(s,0,s.length-1,0);
	}
	/**
	 * sort strings from s[lo] to s[hi] by characters at pos
	 * @param lo starting index
	 * @param hi ending index
	 * @param pos character position
	 */
	private static void sort(String[] s, int lo, int hi, int pos) {
		if (lo >= hi) {
			return;
		}
		int rand = (int) Math.random()*(hi - lo) + lo;
		int i = lo;
		int lt = lo, gt = hi;
		String pivot = s[rand];
		while (i <= hi) {
			int cmp = charat(s[i], pos) - charat(pivot, pos);
			if (cmp == 0) {
				i++;
			} else if (cmp < 0) {
				exch(s, lt++, i++);
			} else {
				exch(s, gt--, i);
			}
		}				
		//System.out.println(String.join(" ", s));
		//strings from pivot[0] to pivot[1] have same char at pos
		//!!!!!
		//charat(pivot, pos + 1) != -1 is WRONG!!!
		//because we are not sure if other strings between lt and gt
		//are longer or shorter than pivot
		if (charat(pivot, pos) != -1) {
			sort(s, lt, gt, pos + 1);
		}
		//remaining strings may have different char at pos
		sort(s, lo, lt - 1, pos);
		sort(s, gt + 1, hi, pos);		
	}
	/**
	 * do 3-way parition on s[lo] to s[hi]
	 * @param s
	 */
	private static int charat(String s, int i) {
		if (s.length() > i) {
			return(s.charAt(i));
		} else {
			return(-1);
		}
	}
	private static void exch(String[] s, int i, int j) {
		String tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	}

}
