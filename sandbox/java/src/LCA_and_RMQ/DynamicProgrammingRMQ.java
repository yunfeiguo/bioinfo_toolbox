package LCA_and_RMQ;

/**
 * use dynamic programming method to
 * calculate range minimum queries
 * @author guoy28
 *
 */
public class DynamicProgrammingRMQ implements RMQ{
    private int N;
    private int[][] result;
    /**
     * should run in O(n^2) time and space
     * use dynamic programming to pre-calculate
     * all possible queries
     * o[i][j] = a[j] < a[o[i][j-1]] ? j : o[i][j-1]
     */
    public DynamicProgrammingRMQ(int[] a) {
        if (a == null || a.length == 0) {
            throw new IllegalArgumentException();
        }
        N = a.length;
        result = new int[N][N];
        result[0][0] = 0;
        for (int j = 1; j < N; j++) {
            result[0][j] = a[result[0][j-1]] <= a[j] ? result[0][j-1] : j;
        }
        for (int i = 1; i < N; i++) {
            result[i][i] = i; //put here for sake of completeness, not really useful
            for (int j = i + 1; j < N; j++) {
                result[i][j] = a[result[i][j - 1]] <= a[j] ? result[i][j - 1] : j;
            }
        }
    }
    public int min(int i, int j) {
        if (i > j || i < 0 || j >= N) {
            throw new IllegalArgumentException();
        }
        return result[i][j];
    }
}
