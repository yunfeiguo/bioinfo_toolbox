package LCA_and_RMQ;
/*
 * implementation of Michael Bender and Martin Farach-Colton's algorithm
 * should achieve O(nlgn) preprocessing and O(1) query time
 * 
 * simple outline of the algorithm
 * for each element i, create lg(n-i) min entries,
 * each entry keeps the index of min element starting
 * from i up to i + 2^(1...lg(n-i)). For minimum range
 * query i...j, we return the index of min of i..i+2^lg(j-i),
 * j-2^lg(j-i)...j.
 * 
 */
public class SparseTableRMQ implements RMQ{
    private int[][] st;
    private int[] pow;
    private int[] a;
    private int N;
    /**
     * initialize a sparse Table n by lgn
     * fill it with index of mins
     * how to fill it efficiently? We use dynamic
     * programming, first fill the first column
     * denoting index of mins of blocks of size 1,
     * then 2nd column, for blocks of size 2. answer
     * for blocks of size 2 can be easily calculated
     * by comparing two adjacent blocks of size 1.
     * @param a the array to be queried
     */
    public SparseTableRMQ(int[] a) {
        this.a = a;
        N = a.length;
        if (N <= 1) {
            return;
        }
        //we get floor of lgN because
        //2*2^lgN is enough for query from 0 to N-1
        int lgN = (int) (Math.log(N)/Math.log(2));
        //pre-calculate powers of 2 to save time
        pow = pow2(lgN);
        st = new int[N][lgN + 1];
        //initialize first column
        for (int i = 0; i < N; i++) {
            st[i][0] = i;
        }
        for (int j = 1; j <= lgN; j++) {
            for (int i = 0; i < N; i++) {
                //not all rows have lgN columns
                //note that i is inclusive
                if (i + pow[j] - 1 >= N) {
                    continue;
                }
                int leftMin = st[i][j-1];
                int rightMin = st[i + pow[j-1]][j - 1];
                //when there are ties, return smallest index
                st[i][j] = a[leftMin] <= a[rightMin] ? leftMin : rightMin;
            }
        }
    }
    /**
     * calculate powers of 2 from 0 to n
     * return an integer array
     * @param n
     * @return
     */
    private int[] pow2(int n) {
        int[] pow = new int[n + 1];
        pow[0] = 1;
        for (int i = 1; i <= n; i++) {
            pow[i] = pow[i - 1] * 2;
        }
        return pow;
    }
    /**
     * return index of min element from i to j,
     * inclusive
     * return -1 for empty array
     */
    public int min(int i, int j) {
        if (i > j || i < 0 || j >= N) {
            throw new IllegalArgumentException();
        }
        if (N == 0) {
            return -1;
        }
        if (i == j || N == 1) {
            return i;
        }
        //find the step where 2*2^step >= j - i + 1
        int stepSize = (int) (Math.log(j - i + 1)/Math.log(2));
        if (pow[stepSize] == (j - i + 1)) {
            stepSize--;
        }
        /*
         * i------------------j
         * i--------i+step-1
         *     j-step+1-------j
         */
        int leftMin = st[i][stepSize];
        int rightMin = st[j - pow[stepSize] + 1][stepSize];
        
        return a[leftMin] <= a[rightMin] ? leftMin : rightMin;
    }

}
