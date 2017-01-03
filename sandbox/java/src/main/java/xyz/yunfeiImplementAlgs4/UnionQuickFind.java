package yunfeiImplementAlgs4;

/**
 * implements union, find interfaces
 * @author Yunfei
 *
 */
public class UnionQuickFind {
	public static void main(String[] args) {
		UnionQuickFind test = new UnionQuickFind(10);
		test.union(0,1);
		test.union(2,1);
		test.union(3,1);
		test.union(2,3);
		test.union(3,9);
		test.union(5,6);
		test.union(6,7);
		System.out.println(test.find(0, 9) == true);
		System.out.println(test.find(0, 7) == false);
		
	}
	private int[] union;
	/**
	 * initilize an array, storing the group id
	 * @param n
	 */
	public UnionQuickFind(int n) {
		union = new int[n];
		for(int i = 0; i < n; i++) union[i] = i;
	}
	/*
	 * connect p and q
	 */
	public void union(int p, int q) {
		int pID = union[p];
		int qID = union[q];
		for (int i = 0; i < union.length; i++) {
			if(union[i] == pID) {
				union[i] = qID;
			}
		}
	}
	public boolean find(int p, int q) {
		return(union[p] == union[q]);
	}
}
