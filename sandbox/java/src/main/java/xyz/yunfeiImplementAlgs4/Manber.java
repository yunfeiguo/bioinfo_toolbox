package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/16/17.
 */

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

import javax.management.RuntimeErrorException;
import java.util.Comparator;

/******************************************************************************
 *  Compilation:  javac Manber.java
 *  Execution:    java Manber < text.txt
 *  Dependencies: StdIn.java
 *
 *  Reads a text corpus from stdin and sorts the suffixes
 *  in subquadratic time using a variant of Manber's algorithm.
 *
 * reference: U Manbers, G Myers, 1989
 ******************************************************************************/
public class Manber implements SuffixArrayInteface{
  private static final int R = 256; //vocabulary size
  private int[] suffixIndices; //sorted suffix index -> original suffix index
  private boolean[] isNewPartition;
  private String s;
  private int n; //length of string
  /**
   * initialize {@code Manber} object
   * @param s
   * @throws NullPointerException if input is null
   */
  public Manber(String s) {
    if (s == null) throw new NullPointerException();
    s = s + '\0'; //assume no '\0' in s
    this.s = s;
    this.n = s.length();
    this.suffixIndices = new int[n];
    //sort array based on 1st char, then each stage, double length of comparison
    int[] counts = new int[R + 2];
    boolean[] isNewPartitionRightHalf = new boolean[n];
    //count
    for (int i = 0; i < n; i++) {
      counts[getChar(s, i) + 2]++;
    }
    //calc boundaries
    for (int i = 1; i < counts.length; i++) {
      counts[i] += counts[i - 1];
    }
    for (int i = 0; i < counts.length; i++) {
      if (counts[i] >= suffixIndices.length)
        break;
      isNewPartitionRightHalf[counts[i]] = true;
    }
    //placement
    for (int i = 0; i < n; i++) {
      int sortedIndex = counts[getChar(s, i) + 1]++;
      suffixIndices[sortedIndex] = i;
    }
    counts = new int[n]; //number of partitions could be N
    int[] reverseSuffixIndices = new int[n]; //original suffix index -> sorted suffix index
    boolean[] isNewPartitionTwoHalf = new boolean[n];
    int[] partitionIndices = new int[n];
    //sort suffixes based on first i chars
    for (int i = 2; i < n*2; i *= 2) {
    //for (int i = 2; i <= 4; i *= 2) {
      sort(suffixIndices, reverseSuffixIndices, counts, partitionIndices, isNewPartitionRightHalf, isNewPartitionTwoHalf, i);
    }
  }

  /**
   * sort suffixes of s based on first k characters assuming
   * suffixes are already sorted based on first k/2 chars
   * @param indices
   * @param reverseSuffixIndices
   * @param counts
   * @param isNewPartitionRightHalf
   * @param isNewPartitionTwoHalf
   * @param k
   */
  private void sort(int[] indices, int[] reverseSuffixIndices,
                    int[] counts, int[] partitionIndices, boolean[] isNewPartitionRightHalf, boolean[] isNewPartitionTwoHalf, int k) {
    if (k < 2) throw new IllegalArgumentException("must have been sorted based on first char");
    int half = k / 2;

    int partitionStartIndex = -1;
    //for each suffix i, we do not need to know its exact location
    //we just need to know where its partition begins
    for (int pos = 0; pos < n; pos++) {
      if (isNewPartitionRightHalf[pos]) {
        partitionStartIndex = pos;
      }
      partitionIndices[indices[pos]] = partitionStartIndex;
      counts[pos] = 0;
    }

    for (int pos = 0; pos < n; pos++) {
      int rightHalfSuffixIndex = indices[pos];
      int twoHalfSuffixIndex = rightHalfSuffixIndex - half;
      if (twoHalfSuffixIndex < 0)
        continue;
      int twoHalfSuffixPartitionStart = partitionIndices[twoHalfSuffixIndex];
      counts[twoHalfSuffixPartitionStart]++;

      int updatedTwoHalfSuffixPosition= twoHalfSuffixPartitionStart + counts[twoHalfSuffixPartitionStart] - 1;
      reverseSuffixIndices[twoHalfSuffixIndex] = updatedTwoHalfSuffixPosition;
      /*
      first set new parition to true
      then reset those elements that are
      not on top of each partition to false
       */
      isNewPartitionTwoHalf[updatedTwoHalfSuffixPosition] = true;
    }
    /*
    within each partition formed by 2 halves
     */
    for (int i = 0; i < n; i++) {
      //if a position is beginning of a new partition, it will
      //always be a beginning of a new partition regardless of
      //k
      indices[reverseSuffixIndices[i]] = i;
    }
    int minUpdatedTwoHalfSuffixPosition = -1;
    int maxUpdatedTwoHalfSuffixPosition = -1;
    int currentPartitionIndex = -1;
    for (int pos = 0; pos < n; pos++) {
      int twoHalfSuffixIndex = indices[pos];
      int rightHalfSuffixIndex = twoHalfSuffixIndex + half;
      if (rightHalfSuffixIndex >= n) {
        maxUpdatedTwoHalfSuffixPosition = minUpdatedTwoHalfSuffixPosition = pos;
        continue;
      }
      if (minUpdatedTwoHalfSuffixPosition == -1) {
        maxUpdatedTwoHalfSuffixPosition = minUpdatedTwoHalfSuffixPosition = pos;
        currentPartitionIndex = partitionIndices[rightHalfSuffixIndex];
      }
      if(isNewPartitionRightHalf[pos] ||
              (rightHalfSuffixIndex < n && currentPartitionIndex != partitionIndices[rightHalfSuffixIndex])) {
        for (int reset = minUpdatedTwoHalfSuffixPosition + 1; reset <= Math.min(n-1, maxUpdatedTwoHalfSuffixPosition); reset++)
          isNewPartitionTwoHalf[reset] = false;
        minUpdatedTwoHalfSuffixPosition = pos;
        maxUpdatedTwoHalfSuffixPosition = minUpdatedTwoHalfSuffixPosition;
        currentPartitionIndex = partitionIndices[rightHalfSuffixIndex];
      } else {
        maxUpdatedTwoHalfSuffixPosition++;
      }
    }
    for (int reset = minUpdatedTwoHalfSuffixPosition + 1; reset <= Math.min(n - 1, maxUpdatedTwoHalfSuffixPosition); reset++)
      isNewPartitionTwoHalf[reset] = false;
    /*
    int minUpdatedTwoHalfSuffixPosition = -1;
    int maxUpdatedTwoHalfSuffixPosition = -1;
    for (int pos = 0; pos < n; pos++) {
      int twoHalfSuffixIndex = indices[pos] - half;
      if (twoHalfSuffixIndex < 0) {
        continue;
      }
      if (minUpdatedTwoHalfSuffixPosition == -1) {
        maxUpdatedTwoHalfSuffixPosition = minUpdatedTwoHalfSuffixPosition = reverseSuffixIndices[twoHalfSuffixIndex];
      }
      if(isNewPartitionRightHalf[pos] ||
              isNewPartitionRightHalf[reverseSuffixIndices[twoHalfSuffixIndex]]) {
        for (int reset = minUpdatedTwoHalfSuffixPosition + 1; reset < Math.min(n, maxUpdatedTwoHalfSuffixPosition); reset++)
          isNewPartitionTwoHalf[reset] = false;
        minUpdatedTwoHalfSuffixPosition = reverseSuffixIndices[twoHalfSuffixIndex];
        maxUpdatedTwoHalfSuffixPosition = minUpdatedTwoHalfSuffixPosition;
      } else {
        maxUpdatedTwoHalfSuffixPosition++;
      }
    }
    for (int reset = minUpdatedTwoHalfSuffixPosition + 1; reset < Math.min(n, maxUpdatedTwoHalfSuffixPosition); reset++)
      isNewPartitionTwoHalf[reset] = false;
    */

    /*
    for (int pos = 0; pos < n; ) {
      int twoHalfSuffixIndex = indices[pos] - half;
      if (twoHalfSuffixIndex < 0) {
        pos++;
        continue;
      }
      int minUpdatedTwoHalfSuffixPosition = reverseSuffixIndices[twoHalfSuffixIndex];
      int maxUpdatedTwoHalfSuffixPosition = minUpdatedTwoHalfSuffixPosition;
      pos++;
      while(pos < n && !isNewPartitionRightHalf[pos]) {
        twoHalfSuffixIndex = indices[pos] - half;
        if (twoHalfSuffixIndex < 0) {
          pos++;
          continue;
        }
        if(isNewPartitionRightHalf[reverseSuffixIndices[twoHalfSuffixIndex]])
          break;
        maxUpdatedTwoHalfSuffixPosition++;
        pos++;
      }
      for (int reset = minUpdatedTwoHalfSuffixPosition + 1; reset <= maxUpdatedTwoHalfSuffixPosition; reset++)
        isNewPartitionTwoHalf[reset] = false;
    }
    */

    for (int i = 0; i < n; i++) {
      //if a position is beginning of a new partition, it will
      //always be a beginning of a new partition regardless of
      //k
      isNewPartitionRightHalf[i] = isNewPartitionRightHalf[i] || isNewPartitionTwoHalf[i];
    }
  }


  /**
   *
   * @param s
   * @param i
   * @return ith char if exists, otherwise -1
   */
  private int getChar(String s, int i) {
    if (i >= 0 && i < s.length())
      return s.charAt(i);
    return -1;
  }

  /**
   *
   * @param i
   * @return original index of ith suffix
   */
  @Override
  public int index(int i) {
    i = i + 1;
    if (i < 0 || i >= suffixIndices.length) throw new IndexOutOfBoundsException();
    //why i + 1 because first suffix is '\0'
    return suffixIndices[i];
  }

  /**
   * @param q
   * @return number of suffixes exactly smaller than q
   */
  @Override
  public int rank(String q) {
    int lo = 1;
    //why lo = 1 because first suffix is '\0'
    int hi = suffixIndices.length - 1;
    while (lo <= hi) {
      int mid = lo + (hi - lo)/2;
      int cmp = compare(s, mid, q);
      if (cmp < 0) {
        lo = mid + 1;
      } else if (cmp > 0) {
        hi = mid - 1;
      } else {
        return mid;
      }
    }
    return lo;
  }

  /**
   * compare s.substring(mid) with q
   * @param s
   * @param i
   * @param q
   */
  private int compare(String s, int i, String q) {
    int l = Math.min(s.length() - index(i), q.length());
    for (int j = 0; j < l; j++) {
      if (s.charAt(j + index(i)) > q.charAt(j)) {
        return 1;
      } else if (s.charAt(j + index(i)) < q.charAt(j)) {
        return -1;
      }
    }
    return s.length() - index(i) - q.length();
  }

  /**
   *
   * @param i
   * @return ith smallest suffix
   */
  @Override
  public String select(int i) {
    if (i < 0 || i >= suffixIndices.length) throw new IndexOutOfBoundsException();
    return s.substring(index(i));
  }
  /**
   * @param i
   * @return length of longest common prefix between i and i-1 suffixes
   */
  @Override
  public int lcp(int i) {
    int l = Math.min(suffixIndices.length - index(i - 1), suffixIndices.length - index(i));
    int j = 0;
    for (; j < l; j++) {
      if (getChar(s, index(i - 1) + j) != getChar(s, index(i) + j)) {
        break;
      }
    }
    return j;
  }

  @Override
  public int length() {
    return s.length() - 1;
  }
  /**
   * Unit tests the {@code SuffixArrayx} data type.
   *
   * @param args the command-line arguments
   */
  public static void main(String[] args) {
    In in = new In(args[0]);
    //String s = in.readAll().replaceAll("\\s+", " ").trim();
    String s = in.readAll().replaceAll("\\W+", " ").trim();
    Manber suffix1 = new Manber(s);
    SuffixArray suffix2 = new SuffixArray(s);
    boolean check = true;
    for (int i = 0; check && i < s.length(); i++) {
      if (suffix1.index(i) != suffix2.index(i)) {
        StdOut.println("suffix1(" + i + ") = " + suffix1.index(i));
        StdOut.println("suffix2(" + i + ") = " + suffix2.index(i));
        String ith = "\"" + s.substring(suffix1.index(i), Math.min(suffix1.index(i) + 50, s.length())) + "\"";
        String jth = "\"" + s.substring(suffix2.index(i), Math.min(suffix2.index(i) + 50, s.length())) + "\"";
        StdOut.println(ith);
        StdOut.println(jth);
        check = false;
      }
    }
    if(!check)
      throw new RuntimeException();
    StdOut.println("  i ind lcp rnk  select");
    StdOut.println("---------------------------");

    for (int i = 0; i < s.length(); i++) {
      int index = suffix1.index(i);
      String ith = "\"" + s.substring(index, Math.min(index + 50, s.length())) + "\"";
      int rank = suffix1.rank(s.substring(index));
      assert s.substring(index).equals(suffix1.select(i));
      if (i == 0) {
        StdOut.printf("%3d %3d %3s %3d  %s\n", i, index, "-", rank, ith);
      }
      else {
        // int lcp  = suffix.lcp(suffix2.index(i), suffix2.index(i-1));
        int lcp  = suffix1.lcp(i);
        StdOut.printf("%3d %3d %3d %3d  %s\n", i, index, lcp, rank, ith);
      }
    }
  }
  /**
   *
   * @param a
   * @param b
   * @return a^b
   */
  private int pow(int a, int b) {
    if (b < 0) throw new IllegalArgumentException("nonnegative must");
    int result = 1;
    for (; b > 0; b--) {
      result *= a;
    }
    return result;
  }
  class IndexComparator implements Comparator<Integer> {
    private int[] a;

    public IndexComparator(int[] a) {
      this.a = a;
    }

    @Override
    public int compare(Integer o1, Integer o2) {
      return a[o1] - a[o2];
    }
  }

  /**
   * use insertion sort to sort suffixes from
   * lo to hi based on k, k+1, k+2, ... characters
   */
  private void insertion(String s, int[] indices, int lo, int hi, int k) {
    for (int i = lo; i <= hi; i++) {
      for (int j = i + 1; j <= hi; j++) {
        int compareIndex = k;
        while(getChar(s, indices[i] + compareIndex) == getChar(s, indices[j] + compareIndex) &&
                indices[i] + compareIndex < s.length() &&
                indices[j] + compareIndex < s.length())
          compareIndex++;
        if (getChar(s, indices[i] + compareIndex) > getChar(s, indices[j] + compareIndex))
          swap(indices, i, j);
      }
    }
  }
  private void swap(int[] a, int i, int j) {
    int tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }
}
