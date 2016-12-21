
public class NumBitsConversion5_5 {
	public static void main(String[] args) {
		//System.out.println(NumBitsConversion5_5.getConvertBits(31,14)==2);
		//System.out.println(NumBitsConversion5_5.getConvertBits(31,31)==0);
		//System.out.println(NumBitsConversion5_5.getConvertBits(-1,1)==31);
		//System.out.println(NumBitsConversion5_5.getConvertBits2(31,14)==2);
		//System.out.println(NumBitsConversion5_5.getConvertBits2(31,31)==0);
		System.out.println(NumBitsConversion5_5.getConvertBits2(-1,1)==31);
	}
	/**
	 * take two numbers, determine # of bits required for conversion
	 * @param a
	 * @param b
	 * @return
	 */
	public static int getConvertBits(int a, int b) {
		int n = a ^ b;
		//count number of 1s in n
		String binary = Integer.toBinaryString(n);
		int i = binary.indexOf('1');
		int count = 0;
		while (i != -1) {
			count++;
			i = binary.indexOf('1', i + 1);
		}
		return(count);
	}
	public static int getConvertBits2(int a, int b) {		
		int count = 0;
		for (int i = a ^ b; i != 0 ; i = i >>> 1) {
			//must use unsigned shift (>>>) here
			System.out.println(i);
			count += i & 1;
		}
		return(count);
	}
}
