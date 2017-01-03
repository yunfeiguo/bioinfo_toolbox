package LCA_and_RMQ;

/**
 * Created by guoy28 on 9/21/16.
 */

import java.util.Arrays;

/**
 * Assume adjacent elements in the input
 * array only differ by +1 or -1, range
 * minimum query (RMQ) problem can be solved
 * in O(n) preprocessing time and O(1) query
 * time following Michael Bender and Martin
 * Colton's algorithm. Specifically, the
 * input array will be divided into blocks
 * of size lgn/2, the total 2n/lgn blocks
 * can be preprocessed and queried in <O(n),
 * O(1)> time by sparse table algorithm.
 * Within each block, RMQ answers can be
 * found quickly by normalizing the block
 * to one of sqrt(n) standard blocks.
 */
public class NormalizedRMQ implements RMQ {
    //array of this size or smaller: use brute force
    private static final int bruteForceCutoff = 3;
    private int[] a;
    private int N;
    private int blockSize;
    //array of intra-block minimums
    private int[] intraBlockMins;
    //array of indeces of intra-block minimums
    private int[] intraBlockMinIdx;
    //RMQ on intraBlockMins
    private RMQ interBlockRMQ;
    /*
     * store all possible RMQ answers for
     * each standardized block. index of
     * a normalized block is encoded by
     * its pattern in binary format. e.g.
     * 011 would be equal to 3, 0101 would
     * be equal to 5.
     */
    private RMQ[] standardizedBlockRMQ;
    //record index of each block in standardizedBlockRMQ
    private int[] block2standardizedBlock;
    //when input array is small, use other implementations of
    //RMQ, when last block is incomplete, use this field
    //to store its answer
    private RMQ smallRMQ;

    public NormalizedRMQ(int[] a) {
        if (a == null || a.length == 0) {
            throw new IllegalArgumentException();
        }
        this.a = a;
        N = a.length;
        /*
         * we don't need to do anything for array of size
         * smaller than 4 because block size will be zero
         */
        if (N <= bruteForceCutoff) {
            smallRMQ = new BruteForceRMQ(a);
            return;
        }
        /*
         * why divided by 2? that will limit number
         * of possible normalized blocks to sqrt(n)
         */
        blockSize = (int) (Math.log(N) / 2);
        //allow last block to be incomplete
        int blockCount = N / blockSize + (N % blockSize == 0 ? 0 : 1);
        intraBlockMinIdx = new int[blockCount];
        intraBlockMins = new int[blockCount];
        for (int i = 0; i < N; i += blockSize) {
            int minIdxCurrentBlock = i;
            for (int j = 1; j < blockSize && (i + j) < N; j++) {
                if (a[i + j] < a[minIdxCurrentBlock]) {
                    minIdxCurrentBlock = i + j;
                }
            }
            intraBlockMinIdx[i / blockSize] = minIdxCurrentBlock;
            intraBlockMins[i / blockSize] = a[minIdxCurrentBlock];
        }
        interBlockRMQ = new SparseTableRMQ(intraBlockMins);
    /*
     * iterate over all blocks (assuming adjacent elements only differ +1
     * or -1), and pre-calculate all possible range min queries for these blocks.
     *
     */
        //we can use RMQ[] to store SparseTableRMQ elements due to inheritance
        standardizedBlockRMQ = new RMQ[1 + (int) (Math.pow(2, blockSize - 1))];
        block2standardizedBlock = new int[blockCount];
        for (int i = 0; i < N; i += blockSize) {
            if (N - i < blockSize) {
                smallRMQ = new SparseTableRMQ(Arrays.copyOfRange(a, i, N));
                //store smallRMQ as last element of standardizedBlockRMQ
                //so we don't have to treat it differently
                block2standardizedBlock[i / blockSize] = standardizedBlockRMQ.length - 1;
                standardizedBlockRMQ[standardizedBlockRMQ.length - 1] = smallRMQ;
                break;
            }
            int standardizedBlockIdx = 0;
            for (int j = 1; j < blockSize; j++) {
                //assume adjacent elements differ by +1 or -1
                standardizedBlockIdx <<= 1;
                if (a[i + j] - a[i + j - 1] == 1) {
                    standardizedBlockIdx++;
                }
            }
            if (standardizedBlockRMQ[standardizedBlockIdx] == null) {
                standardizedBlockRMQ[standardizedBlockIdx] = new SparseTableRMQ(Arrays.copyOfRange(a, i, i + blockSize));
            }
            block2standardizedBlock[i / blockSize] = standardizedBlockIdx;
        }
    }


    /**
     * return index of min index between i and j
     * i must be no larger than j
     *
     * @param i
     * @param j
     * @return
     * @throws IllegalArgumentException
     */
    public int min(int i, int j) {
        if (i < 0 || j >= this.N || i > j) {
            throw new IllegalArgumentException();
        }
        if (this.N <= bruteForceCutoff) {
            return this.smallRMQ.min(i, j);
        }
        /*
        2 cases:
        1) i and j are in the same block
        we use StandardidizedBlockRMQ to solve it
        2) i and j are in different blocks
        min can be found by comparing i->block boundary,block boundary->j
        and any blocks (could be none) in between.
         */
        if (i / blockSize == j / blockSize) {
            //same block
            return intraBlockMin(i, j);
        }
        int leftMin = intraBlockMin(i, (1 + i / blockSize) * blockSize - 1);
        int rightMin = intraBlockMin(j - j % blockSize, j);
        int marginMin = a[leftMin] <= a[rightMin] ? leftMin : rightMin;
        if (j / blockSize - i / blockSize > 1) {
            //there are blocks in between
            int middleMin = intraBlockMinIdx[interBlockRMQ.min(i / blockSize + 1, j / blockSize - 1)];
            if (a[marginMin] != a[middleMin]) {
                return a[marginMin] < a[middleMin] ? marginMin : middleMin;
            } else {
                //when equal, return smaller index
                return a[leftMin] == a[middleMin] ? leftMin : middleMin;
            }
        } else {
            return marginMin;
        }
    }

    /**
     * assume i and j are within same block
     * return index of the min element between
     * i and j
     *
     * @param i
     * @param j
     * @return
     */
    private int intraBlockMin(int i, int j) {
        if (i / blockSize != j / blockSize) {
            throw new IllegalArgumentException("not in the same block!");
        }
        return (i - i % blockSize) +
                standardizedBlockRMQ[block2standardizedBlock[i / blockSize]].min(i % blockSize,
                        j % blockSize);
    }

    /**
     * unit tests
     * test the following cases
     * 1) size of 1 array
     * 2) size of 8 array, block size 1
     * 3) size of 18 array, inter block query (no block in middle), intra block query
     * inter block query (1 block in middle)
     *
     * @param args
     */
    public static void main(String[] args) throws Exception {
        RMQ test = new NormalizedRMQ(new int[]{0});
        if (test.min(0, 0) != 0) {
            throw new Exception();
        }
        test = new NormalizedRMQ(new int[]{0, 1, 2, 1, 0, -1, 0, 1});
        if (test.min(0, 7) != 5) {
            throw new Exception();
        }
        test = new NormalizedRMQ(new int[]{0, 1, 2, 1, 0, -1, 0, 1, 0, 1, 2, 1, 0, -1, 0, 1, 2, 3});
        if (test.min(4, 5) != 5) {
            throw new Exception();
        }
        if (test.min(5, 6) != 5) {
            throw new Exception();
        }
        if (test.min(2, 17) != 5) {
            throw new Exception();
        }
    }
}
