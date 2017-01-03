package yunfeiImplementAlgs4;

public class CalcQuickSortTime {
	public static void main(String[] args) {
		int[] n = {100, 1000, 10000};
		for (int i : n) {
			System.out.println("n:" + i);
			//approximation
			double approx = 2*Math.log(i)*i;
			//exact
			double exact = 0;
			if (i < 2) {
				System.out.println("calc by hand");
			} else {				
				for (int j = 2; j <= i; j++) {
					exact = 2 + (j + 1)/j*exact;
				}
			}
			double ratio = approx/exact;
			System.out.println(ratio);
		}
	}
}
