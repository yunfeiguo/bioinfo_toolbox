package programming_assignments.HW5_BurrowsWheeler.HW5_BurrowsWheeler;

/**
 * Created by guoy28 on 5/7/17.
 */
public class CircularSuffixArray {
  private static int R = 256;
  private int[] indices;
  // circular suffix array of s
  public CircularSuffixArray(String in) {
    if (in == null) {
      throw new NullPointerException();
    }
    int n = in.length();
    /*sort by stages
    first stage, sort by 1st char
    second stage, sort by first 2 chars
    third stage, sort by first 4 chars
    ...
    in total, finish sorting in lgN stages
    */
    char[] a = in.toCharArray();
    int[] indices = new int[n];
    this.indices = indices;
    for (int i = 0; i < n; i++) {
      indices[i] = i;
    }
    int[] boundaries = sort(a, indices);
    /*
    AB CXX
    AB TYY
    AC CTT
    AT ABC
    AT ABT
    AT ATT
    AT ACC
    old boundaries: 0 0 2 3 3 3 3 ...
    new boundaries: 0 1 2 3 3 5 4 ...
     */
    //step size
    for (int s = 2; s <= n; s += s) {
      int[] newBoundaries = new int[n];
      int[] newIndices = new int[n];
      //counts of elements in each bucket marked by s/2 prefix
      int[] leftBoundaryCounts = new int[n];
      //keep track of new boundaries relative to s/2 prefix boundaries
      int[] leftNewRelativeBoundary = new int[n];
      for (int i = 0; i < n; i++) {
        /*
        traverse indices in sequential order
        record start of each bucket
        for each suffix, it determines
        location of -s/2 suffix
         */
        int currentIndex = indices[i];
        int prefixIndex = (n + currentIndex - s/2) % n;
        //boundaries of suffix i
        int currentRightBoundary = boundaries[currentIndex];
        //boundaries of suffix i - s/2
        int currentLeftBoundary = boundaries[prefixIndex];
        int currentBoundaryCount = leftBoundaryCounts[currentLeftBoundary] - leftNewRelativeBoundary[currentLeftBoundary];
        /*
        assume prefixIndex is in the same bucket as previous element with same prefix
        then their suffices should be the same.
         */
        int firstIndex = newIndices[(currentLeftBoundary + leftNewRelativeBoundary[currentLeftBoundary]) % n];
        if (currentBoundaryCount > 0 && boundaries[(firstIndex + s/2) % n] == currentRightBoundary) {
          newBoundaries[prefixIndex] = currentLeftBoundary + leftNewRelativeBoundary[currentLeftBoundary];
        } else {
          //got new bucket
          leftNewRelativeBoundary[currentLeftBoundary] = leftBoundaryCounts[currentLeftBoundary];
          newBoundaries[prefixIndex] = currentLeftBoundary + leftNewRelativeBoundary[currentLeftBoundary];
        }
        currentBoundaryCount = leftBoundaryCounts[currentLeftBoundary] - leftNewRelativeBoundary[currentLeftBoundary];
        leftBoundaryCounts[currentLeftBoundary]++;
        newIndices[(newBoundaries[prefixIndex] + currentBoundaryCount) % n] = prefixIndex;
      }
      System.arraycopy(newBoundaries, 0, boundaries, 0, boundaries.length);
      System.arraycopy(newIndices, 0, indices, 0, indices.length);
    }
  }

  /**
   * sort indices based on array of chars
   * @param a
   * @param indices
   */
  private int[] sort(char[] a, int[] indices) {
    //time: O(n)
    int n = a.length;
    int[] count = new int[R + 1];
    int[] boundaries = new int[n];
    //count each char
    for (char c : a) {
      count[c + 1]++;
    }
    //calc boundaries
    for (int i = 1; i < count.length; i++) {
      count[i] += count[i - 1];
    }

    for (int i = 0; i < n; i++) {
      boundaries[indices[i]] = count[a[i]];
    }
    int[] newIndices = new int[n];
    for (int i = 0; i < n; i++) {
      newIndices[count[a[i]]++] = indices[i];
    }
    System.arraycopy(newIndices, 0, indices, 0, n);
    return boundaries;
  }

   // length of s
  public int length() {
    return indices.length;
  }
              // returns index of ith sorted suffix
  public int index(int i) {
    if (i < 0 || i >= length())
      throw new IndexOutOfBoundsException();
    return indices[i];
  }
  public static void main(String[] args) {
    String s = "ABRACADABRA!";
    CircularSuffixArray a = new CircularSuffixArray(s);
    for (int i = 0; i < a.length(); i++) {
      System.out.print(a.index(i));
      System.out.print(" ");
    }
    System.out.println("");
    System.out.println("expect: ");
    System.out.println("11 10 7 0 3 5 8 1 4 6 9 2");
  }// unit testing of the methods (optional)
}
