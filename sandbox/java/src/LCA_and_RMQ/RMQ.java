package LCA_and_RMQ;

/*
 * interface for range minimum query problem
 * 
 */
public interface RMQ {
    /**
     * returns index of min element in array[i...j]
     * when there are ties, return the first index
     * @param i
     * @param j
     * @return
     */
    public int min(int i, int j);
}
