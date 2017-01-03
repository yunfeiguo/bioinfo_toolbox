package yunfeiImplementAlgs4;

public class Quick3String {
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
		int[] pivot = partition(s,lo,hi,pos);
		if (pivot[0] == -1) {
			return;
		}
		//System.out.println(String.join(" ", s));
		//strings from pivot[0] to pivot[1] have same char at pos
		sort(s, pivot[0], pivot[1], pos + 1);
		//remaining strings may have different char at pos
		sort(s, lo, pivot[0] - 1, pos);
		sort(s, pivot[1] + 1, hi, pos);		
	}
	/**
	 * do 3-way parition on s[lo] to s[hi]
	 * @param s
	 */
	private static int[] partition(String[] s, int lo, int hi, int pos) {
		int rand = (int) Math.random()*(hi - lo) + lo;
		int i = lo;
		String pivot = s[rand];
		boolean notEnd = false;
		while (i <= hi) {
			if (s[i].length() > pos) {
				notEnd = true;
			}
			int cmp = compare(s[i],pivot,pos);
			if (cmp == 0) {
				i++;
			} else if (cmp < 0) {
				exch(s, lo++, i++);
			} else {
				exch(s, hi--, i);
			}
		}		
		if (!notEnd) {
			return(new int[]{-1,0});
		} else {
			return(new int[] {lo,hi});
		}		
	}
	private static int compare(String a, String b, int pos) {
		if (a.length() > pos) {
			if (b.length() > pos) {
				return(a.charAt(pos) - b.charAt(pos));
			} else {
				return(1);
			}
		} else {
			if (b.length() > pos) {
				return(-1);
			} else {
				return(0);
			}
		}
	}
	private static void exch(String[] s, int i, int j) {
		String tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	}

}
