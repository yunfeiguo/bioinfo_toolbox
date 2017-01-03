package LCA_and_RMQ;

/*
 * brute force method for range minimum query problem
 */
public class BruteForceRMQ implements RMQ{
    private int N;
    private int[][] result;
    /**
     * take an array of integers and initialize query object
     * create a n^2 table to store results of all possible queries
     * the table can be generated in O(n^3) time
     * @param a
     */
    public BruteForceRMQ(int[] a) {
        this.N = a.length;
        this.result = new int[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                this.result[i][j] = min(a, i, j);
            }
        }
    }
    /**
     * return index of min element in a[i..j]
     * when there are ties, return first one
     * @param a
     * @param i
     * @param j
     * @return
     */
    private int min(int[] a, int i, int j) {
        if (i > j) {
            throw new IllegalArgumentException();
        }
        int minIndex = i;
        for (int k = i + 1; k <= j; k++) {
            if (a[k] < a[minIndex]) {
                minIndex = k;
            }
        }
        return minIndex;
    }
    /**
     * return index of min element in a[i..j]
     * when there are ties, return the first index
     * @param i
     * @param j
     * @return
     */
    public int min(int i, int j) {
        if (i > j || i < 0 || j >= N) {
            throw new IllegalArgumentException();
        }
        return result[i][j];
    }
}
