package yunfeiImplementAlgs4;

public class QuickInt {
        public int[] quickSort(int[] array) {
                if (array == null) {
                    return array;
                }
                quickSort(array, 0, array.length - 1);
                return array;
            }
            private void quickSort(int[] a, int lo, int hi) {
                //lo: lower boundary (inclusive)
                //hi: upper boundary (inclusive)
                if (lo >= hi) {
                    return;
                }
                int pivot = partition(a, lo, hi);
                quickSort(a, lo, pivot - 1);
                quickSort(a, pivot + 1, hi);
            }
            private int partition(int[] a, int lo, int hi) {
                int pivot = (int) (Math.random() * (hi - lo + 1)) + lo;
                exch(a, pivot, lo);
                int i = lo; //left pointer
                int j = hi + 1; //right pointer
                while (true) {
                    while(a[++i] < a[lo]) {
                        if (i == hi) {
                            break;
                        }
                    }
                    while(a[lo] < a[--j]) {
                        if (j == lo) {
                            break;
                        }
                    }
                    if (i >= j) {
                        break;
                    } else {
                        exch(a, i, j);
                    }
                }   
        exch(a, lo, j); 
                return j;
            }   
            private void exch(int[] a, int i, int j) {
                int tmp = a[i];
                a[i] = a[j];
                a[j] = tmp;
            }        
}
