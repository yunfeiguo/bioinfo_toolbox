package programming_assignments.HW5_BurrowsWheeler.HW5_BurrowsWheeler2;
/* this implementation uses Manber and Myer's method
 * to sort N circular suffix arrays in O(NlgN) time.
 * the speedup is a result of utilizing information
 * gained during sorting prefixes
 */
public class CircularSuffixArray {
    private final String s;  
    private final int alphabetSize = 256;
    private final int N;
    private int[] array;
    /**
     * sort N circular suffixes in linear space
     * and O(NlgN) time
     * @param s
     */
    public CircularSuffixArray(String s) {
        // circular suffix array of s
        if (s == null) {
            throw new java.lang.NullPointerException();
        }
        this.s = s;
        this.N = s.length();
        this.array = new int[s.length()];
        for (int i = 0; i < s.length(); i++) {
            this.array[i] = i;
        }
        sort();
    }
    public int length() {
        // length of s
        return s.length();
    }
    public int index(int i) {
        // returns index of ith sorted suffix
        if (i < 0 || i >= length()) {
            throw new java.lang.IndexOutOfBoundsException();
        }
        return this.array[i];
    }
    /***************************************************/
    /* perform sorting in lg(N+1) stages,
     * first stage sorts CSA by first char,
     * second stage sorts CSA by first 2 chars,
     * 3rd stage sorts CSA by first 4 chars,
     * ...
     * 
     * utilize the fact that ith suffix's first H
     * chars are same as (i-H)th suffix's 2nd H chars
     * 
     * first round of sorting is done by radix sort
     */
    private void sort() {        
        //count[i] records how many elements in
        //bucket rank[i] have been moved
        int[] count = new int[N];
        //bucket number        
        int[] rank = radixSort();
        //inverse of this.array
        int[] position = new int[N];        
        //array has been sorted based on first char
        //by radixSort()     
        boolean[] isFirstInHBucket = new boolean[N];
        for (int i = 0; i < N; i++) {
            position[array[i]] = i;
        }
        for (int i = 0; i < N; i++) {
            isFirstInHBucket[i] = position[i] == rank[i]; 
        }
        for (int h = 1; h < N; h += h) {
            //initialize element position to be bucket position
            for (int i = 0; i < N; i++) {
                position[i] = rank[i];
                count[i] = 0;
            }
            //update array[i]'s position based on
            //array[i] + h's bucket
            //recall that we are visiting the H+1 to 2H characters            
            for (int i = 0; i < N; i++) {
                int leftShiftHIndex = array[i] - h;
                leftShiftHIndex = leftShiftHIndex < 0? leftShiftHIndex + N : leftShiftHIndex;
                if (rank[leftShiftHIndex] == -1 || leftShiftHIndex == -1) {
                    System.out.println(1);
                }
                position[leftShiftHIndex] += count[rank[leftShiftHIndex]]++;
            }
            //update array and isFirstInHBucket
            for (int i = 0; i < N; i++) {
                array[position[i]] = i;                
            }   
            //update isFirstInHBucket 
            //originalIndex + h's rank marks the H+1 to 2H bucket
            //every time it changes, we mark isFirstInHBucket true
            int originalIndex = 0;
            int mostRecentHRank = -1; //make sure it is different from first rank
            int mostRecent2HRank = -1; //it's possible that first H chars are different, but second H chars are the same, they should still belong to different buckets
            int currentHRank = -1;
            int current2HRank = -1;
            for (int i = 0; i < N; i++) {
                originalIndex = array[i];
                current2HRank = rank[(originalIndex + h) % N];
                currentHRank = rank[originalIndex];
                if (current2HRank != mostRecent2HRank || currentHRank != mostRecentHRank) {
                    isFirstInHBucket[originalIndex] = true;
                } else {
                    isFirstInHBucket[originalIndex] = false;
                }
                mostRecentHRank = currentHRank;
                mostRecent2HRank = current2HRank;
            }
            //update rank i+h based on rank of i
            //we need to update rank based on sorted
            //order, otherwise it's hard to track
            //most recent rank
            int mostRecentRank = -1;
            for (int i = 0; i < N; i++) {
                originalIndex = array[i];                               
                if (isFirstInHBucket[originalIndex]) {
                    mostRecentRank = position[originalIndex];
                }
                rank[originalIndex] = mostRecentRank;
            }            
        }
    }
    /*
     * sort CSA by first char using radix sort
     */
    private int[] radixSort() {
        int[] count = new int[alphabetSize + 1];
        int[] rank = new int[N];
        //count[i] records count of char(i-1)
        for (int i = 0; i < N; i++) {
            count[s.charAt(i) + 1]++;
        }
        //count[i] records start position of i+1 element
        for (int i = 1; i <= alphabetSize; i++) {
            count[i] += count[i - 1];
        }
        //get rank for each csa
        for (int i = 0; i < N; i++) {
            rank[i] = count[s.charAt(i)];
        }
        //convert count to actual indices
        for (int i = 0; i < N; i++) {
            array[count[s.charAt(i)]++] = i;            
        }   
        return rank;
    }
    /***************************************************/
    public static void main(String[] args) {
        // unit testing of the methods (optional)
        CircularSuffixArray test = new CircularSuffixArray("ABRACADABRA!");
        //check some of the sorted indices
        System.out.println(test.index(0) == 11 && test.index(1) == 10 && test.index(11) == 2);
        CircularSuffixArray test2 = new CircularSuffixArray("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
        System.out.println(test2.index(0));
    }
}

