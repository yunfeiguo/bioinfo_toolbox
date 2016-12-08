public class Solution {
    /**
     * @param A, B: Two string.
     * @return: the length of the longest common substring.
     */
    public int longestCommonSubstring(String A, String B) {
        // write your code here
        if(A == null || B == null || A.length() == 0 || B.length() == 0) {
            return(0);
        }
        //use DP
        //lcs[i,j] store the length of lcs ending at i-1 in A, j-1 in B
        //lcs[i,j] = lcs[i-1,j-1] if A[i-1] == B[j-1]
        //lcs[i,j] = 0 otherwise
        //lcs must end somewhere in A, and B, so we search for the max
        //and return it
        int ret = 0;
        int nA = A.length();
        int nB = B.length();
        int[][] lcs = new int[nA + 1][nB + 1];
        for (int i = 1; i <= nA; i++) {
            for (int j = 1; j <= nB; j++) {
                if(A.charAt(i - 1) == B.charAt(j - 1)) {
                    lcs[i][j] = lcs[i - 1][j - 1] + 1;
                    ret = Math.max(lcs[i][j], ret);
                }
            }
        }
        return(ret);
    }
}
