
public class SwapBit5_6 {
	public static void main(String[] args) {
		System.out.println(SwapBit5_6.swap(1)==2);
		System.out.println(SwapBit5_6.swap(0)==0);
		System.out.println(SwapBit5_6.swap(-1)==-1);
		System.out.println(SwapBit5_6.swap2(1)==2);
		System.out.println(SwapBit5_6.swap2(0)==0);
		System.out.println(SwapBit5_6.swap2(-1)==-1);
	}
	public static int swap(int n) {
		int len = 32;//assume 32 bit
		for(int i = 1; i < len; i+=2) {
			boolean cur = getBit(n,i);
			boolean prev = getBit(n, i - 1);
			n = setBit(n, i, prev);
			n = setBit(n, i - 1, cur);
			//System.out.println(n);
		}
		return(n);
	}
	public static int swap2(int n) {
		//this only costs five instructions!
		return( ((n & 0b10101010101010101010101010101010) >>> 1) |
				((n & 0b01010101010101010101010101010101) << 1) );
	}
	public static boolean getBit(int n, int i) {
		return((n & (1 << i)) != 0);
	}
	public static int setBit(int n, int i, boolean t) {
		if(t) {
			return(n | (1 << i));
		} else {
			int mask = ~(1 << i);
			return(n & mask);
		}
	}

}
