package yunfeiImplementAlgs4;

public class Sort2distinct {
	public static void main(String[] args) {
		Integer[] a = new Integer[4];
		a[0] = 1;
		a[1] = 2;
		a[2] = 1;
		a[3] = 2;
		/*a[4] = 2;
		a[5] = 1;*/
		Sort2distinct.sort(a);
		for(int i = 0; i < a.length; i++) {
			System.out.println(a[i]);
		}
	}
	public static void sort(Comparable[] a) {
		int lo = 0, hi = a.length - 1;
		int i = 0;		
		while (i <= hi) {
			int cmp = a[i].compareTo(a[hi]);
			if ( cmp > 0) {
				exch(a, i, hi--);								
			} else if (cmp < 0) {
				exch(a, i++, lo++);				
			} else i++;			
		}		
	}
	private static void exch(Comparable[] a, int i, int j) {
		Comparable tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	}

}
