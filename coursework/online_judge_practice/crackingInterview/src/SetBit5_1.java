
public class SetBit5_1 {
	public static void main(String[] args) {
		int n = Integer.parseInt("10000000000",2);
		int m = Integer.parseInt("10101", 2);
		int ret = Integer.parseInt("10001010100",2);
		System.out.println(ret == SetBit5_1.setBit(n, m, 2, 6));
	}
	/**
	 * set bits from i to j in n to m
	 * @param n
	 * @param m
	 * @param i
	 * @param j
	 * @return
	 */
	public static int setBit(int n, int m, int i, int j) {
		m = m << i;
		int mask = (1 << (j - i + 1)) - 1;
		mask = mask << i;
		n = n & (~ mask);
		n = n | m;
		return(n);
	}

}
