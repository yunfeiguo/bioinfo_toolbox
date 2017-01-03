package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.Stopwatch;

public class TestQuick {
        public static void main(String[] args) {
            QuickCoco test = new QuickCoco();
            QuickInt testYG = new QuickInt();
            int N = 20000000;
            int[] a = new int[N];
            Integer[] b = new Integer[N];
            int[] c = new int[N];
            for (int i = 0; i < N; i++) {
                a[i] = N - i;
                b[i] = N - i;
                c[i] = N - i;
            }             
            Stopwatch timer = new Stopwatch();
            Quick.sort(b);
            System.out.println("Yunfei time: " + timer.elapsedTime());
            for (int i = 0; i < 5; i++) {
                System.out.println(b[i]);
            }
            timer = new Stopwatch();
            c = testYG.quickSort(c);
            System.out.println("Yunfei int time: " + timer.elapsedTime());
            for (int i = 0; i < 5; i++) {
                System.out.println(c[i]);
            }
            timer = new Stopwatch();
            a = test.quickSort(a);
            System.out.println("Coco time: " + timer.elapsedTime());
            for (int i = 0; i < 5; i++) {
                System.out.println(a[i]);
            }
        }
}
