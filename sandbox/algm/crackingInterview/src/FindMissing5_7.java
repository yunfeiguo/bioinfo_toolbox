
public class FindMissing5_7 {
	public static void main(String[] args) {		
		System.out.println(FindMissing5_7.findMissing(new int[]{0,1,2}) == 3);
		System.out.println(FindMissing5_7.findMissing(new int[]{0,1,2,4}) == 3);
		System.out.println(FindMissing5_7.findMissing(new int[]{0,1,4,2,5}) == 3);
		System.out.println(FindMissing5_7.findMissing(new int[]{6,0,5,4,3,2}) == 1);
	}
	public static int findMissing(int[] a) {
		int ret = 0;
		int n = a.length;
		if (n == 0) return(-1);
		for(int i = 0; i < 32; i++) {
			//assume 32bit
			for(int j = 0; j < a.length; j++) {
				ret ^= a[j] & (1 << i);
			}
		}
		if( n % 2 == 0) {
			ret ^= n + 1;
			if ( ((n + 2 - 2)/2) % 2 == 0) {
				ret ^= 1;//flip last bit
			}
		} else {
			if ( ((n + 1 - 2)/2) % 2 == 0) {
				ret ^= 1;
			}
		}
		return(ret);
	}

}
