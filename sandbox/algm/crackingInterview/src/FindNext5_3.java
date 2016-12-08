
public class FindNext5_3 {
	public static void main(String[] args) {
		int x = FindNext5_3.findNext(7);
		int y = FindNext5_3.findNext(6);
		
		System.out.println( x == 11);		
		System.out.println( y == 9);
	}
	//find Next smallest is ommitted here
	public static int findNext(int n) {
		if(n < 0) return(-1);
		int zeroIdx = 0;
		int countOne = 0;
		//find first one
		while(!getBit(n, zeroIdx)) zeroIdx++;
		//find first zero after one
		while(getBit(n, zeroIdx)) {
			zeroIdx++;
			countOne++;
		}
		//0 => 1
		n = setBit(n, zeroIdx--, true);		
		n = setBit(n, zeroIdx--, false);
		countOne--;
		for(int i = zeroIdx; i >= countOne; i--) {
			n = setBit(n, i, false);
		}
		for(int i = 0; i < countOne; i++) {
			n = setBit(n, i, true);
		}
		return(n);		
	}
	public static boolean getBit(int n, int i) {
		return((n & (1 << i)) > 0);
	}
	public static int setBit(int n, int i, boolean flag) {
		if(flag) {
			return(n | (1 << i));
		} else {
			int mask = 1 << i;
			return(n & (~mask));
		}
	}

}
