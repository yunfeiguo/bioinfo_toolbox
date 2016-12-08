public class StringCompare {
    public int minimumTotal(int[][] triangle){
        if(triangle == null || triangle.length == 0){
            return Integer.MAX_VALUE;
        }

        int[][] result = new int[triangle.length][triangle[0].length];
        int globalMin = Integer.MAX_VALUE;
        for(int i = 0; i < triangle.length; i++){
            for(int j = 0; j <= i; j++){
                // initialization
                if(i == 0){
                    result[i][j] = triangle[i][j];
                }else if(j == 0){
                    result[i][j] = triangle[i][j] + result[i - 1][j];
                }else if(j == i){
                    result[i][j] = triangle[i][j] + result[i - 1][j - 1];
                }else{
                    // dp
                    result[i][j] = Math.min(result[i - 1][j - 1], result[i - 1][j]) + triangle[i][j];
                }
                if(i == triangle.length - 1){
                    globalMin = Math.min(globalMin, result[i][j]);
                }
            }
        }
        return globalMin;
    }
    public static void main(String[] args) {
        int[][] triangle = new int[2][];
        triangle[0] = 1;
        triangle[1][0] = 2;
        triangle[1][1] = 3;
        StringCompare test = new StringCompare();
        System.out.println(test.minimumTotal(triangle));
    }
}
