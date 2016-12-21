package LCA_and_RMQ;
import java.util.Random;
import edu.princeton.cs.algs4.Stopwatch;

public class TestRMQ {
    public static void main(String[] args) throws Exception {
        Random random = new Random();
        random.setSeed(11);
        Stopwatch watch = new Stopwatch();
        int N = 1000;
        int[] a = new int[N];
        //randomly generate an array of integers
        for (int i = 0; i < N; i++) {
            a[i] = random.nextInt(N);
        }
        double bfPreprocessTime = 0;
        double dpPreprocessTime = 0;
        double stPreprocessTime = 0;
        bfPreprocessTime -= watch.elapsedTime();
        RMQ bf = new BruteForceRMQ(a);
        bfPreprocessTime += watch.elapsedTime();
        dpPreprocessTime -= watch.elapsedTime();
        RMQ dp = new DynamicProgrammingRMQ(a);
        dpPreprocessTime += watch.elapsedTime();
        stPreprocessTime -= watch.elapsedTime();
        RMQ st = new SparseTableRMQ(a);
        stPreprocessTime += watch.elapsedTime();

        double bfQueryTime = 0;
        double dpQueryTime = 0;
        double stQueryTime = 0;
        //randomly test some ranges
        int nTest = 1000;
        int count = nTest;
        while(count > 0) {
            int i = random.nextInt(N);
            int j = i + random.nextInt(N - i);
            bfQueryTime -= watch.elapsedTime();
            int answer = bf.min(i,j);
            bfQueryTime += watch.elapsedTime();

            System.out.println("testing " + i + " " + j);
            dpQueryTime -= watch.elapsedTime();
            if(answer != dp.min(i,j)) {
                throw new Exception("For dp: Expected " + answer + " Got " + dp.min(i, j));
            }
            dpQueryTime += watch.elapsedTime();
            stQueryTime -= watch.elapsedTime();                        
            if(answer != st.min(i,j)) {
                throw new Exception("For st: Expected " + answer + " Got " + st.min(i, j));
            }
            stQueryTime += watch.elapsedTime();
            count--;
        }
        System.out.println("array size: " + N);
        System.out.println(nTest + " tests done");
        

        System.out.println("Brute force preprocessing time (s): " + bfPreprocessTime);
        System.out.println("Dynamic programming preprocessing time (s): " + dpPreprocessTime);
        System.out.println("Sparse table preprocessing time (s): " + stPreprocessTime);
        
        System.out.println("Brute force query time (s): " + bfQueryTime);
        System.out.println("Dynamic programming query time (s): " + dpQueryTime);
        System.out.println("Sparse table query time (s): " + stQueryTime);
        
    }
}
