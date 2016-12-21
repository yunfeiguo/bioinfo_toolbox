package yunfeiImplementAlgs4;

public class QuickCoco {
    public int[] quickSort(int[] array) {
        if(array == null || array.length <= 1){
            return array;
        }
        quickSortHelper(array, 0, array.length - 1);
        return array;
    }
    private void quickSortHelper(int[] array, int left, int right){
        // termination
        if(left >= right){
            return;
        }
        int pivotIndex = partition(array, left, right);
        quickSortHelper(array, left, pivotIndex - 1);
        quickSortHelper(array, pivotIndex + 1, right);
        return;
    }
    private int partition(int[] array, int left, int right){
        int pivotIndex = getPivot(left, right);
        int pivotElement = array[pivotIndex];
        // partition
        // swap pivot to the end
        swap(array, pivotIndex, right);
        // from left to right - 1, swap everything larger than pivot to the right
        int leftIndex = left;
        int rightIndex = right - 1;
        while(leftIndex <= rightIndex){
            if(array[leftIndex] < pivotElement){
                leftIndex++;
            }else if(array[rightIndex] > pivotElement){
                rightIndex--;
            }else{
                swap(array, leftIndex, rightIndex);
                leftIndex++;
                rightIndex--;
            }
        }
        // swap back pivot
        swap(array, leftIndex, right);
        // return pivotIndex
        return leftIndex;
    }
    private int getPivot(int left, int right){
        return left + (int) (Math.random() * (right - left + 1));
    }
    private void swap(int[] array, int left, int right){
        int tmp = array[left];
        array[left] = array[right];
        array[right] = tmp;
    }

}
