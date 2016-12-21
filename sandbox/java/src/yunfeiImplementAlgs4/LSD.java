package yunfeiImplementAlgs4;
public class LSD {
	public static void sort(String[] a, int W) {
		int N = a.length;
		String[] aux = new String[N];
		int R = 256;
		
		for (int i = W - 1; i >= 0; i--) {
			int[] count = new int[R+1];
			//count occurrences
			for (int j = 0; j < N; j++) {
				count[a[j].charAt(i) + 1]++;
			}
			//calculate starting index
			for (int j = 0; j < R; j++){
				count[j + 1] += count[j];
			}
			//distribute results
			//does not change relative position
			//of strings with same character at i
			for (int j = 0; j < N; j++) {
				aux[count[a[j].charAt(i)]++] = a[j];
			}
			//copy back
			for (int j = 0; j < N; j++) {
				a[j] = aux[j];
			}
		}
	}
}
