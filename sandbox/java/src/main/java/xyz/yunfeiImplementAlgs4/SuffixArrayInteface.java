package yunfeiImplementAlgs4;

/**
 * Created by guoy28 on 3/19/17.
 */
public interface SuffixArrayInteface {
    /**
     * Returns the length of the input string.
     *
     * @return the length of the input string
     */
    public int length();

    /**
     * Returns the index into the original string of the <em>i</em>th smallest suffix.
     * That is, {@code text.substring(sa.index(i))} is the <em>i</em>th smallest suffix.
     *
     * @param i an integer between 0 and <em>n</em>-1
     * @return the index into the original string of the <em>i</em>th smallest suffix
     * @throws java.lang.IndexOutOfBoundsException unless {@code 0 <= i < n}
     */
    public int index(int i);

    /**
     * Returns the length of the longest common prefix of the <em>i</em>th
     * smallest suffix and the <em>i</em>-1st smallest suffix.
     *
     * @param i an integer between 1 and <em>n</em>-1
     * @return the length of the longest common prefix of the <em>i</em>th
     * smallest suffix and the <em>i</em>-1st smallest suffix.
     * @throws java.lang.IndexOutOfBoundsException unless {@code 1 <= i < n}
     */
    public int lcp(int i);

    /**
     * Returns the <em>i</em>th smallest suffix as a string.
     *
     * @param i the index
     * @return the <em>i</em> smallest suffix as a string
     * @throws java.lang.IndexOutOfBoundsException unless {@code 0 <= i < n}
     */
    public String select(int i);

    /**
     * Returns the number of suffixes strictly less than the {@code query} string.
     * We note that {@code rank(select(i))} equals {@code i} for each {@code i}
     * between 0 and <em>n</em>-1.
     *
     * @param query the query string
     * @return the number of suffixes strictly less than {@code query}
     */
    public int rank(String query);
}
