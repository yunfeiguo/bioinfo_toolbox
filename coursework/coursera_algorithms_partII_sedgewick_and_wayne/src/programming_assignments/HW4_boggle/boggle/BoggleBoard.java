package programming_assignments.HW4_boggle.boggle;

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdRandom;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * Created by guoy28 on 4/5/17.
 */
public class BoggleBoard {
  private char[][] board;
  private static final char[] hasbroDice = new char[]{'A','T','E','E','A','P','Y','O','T','I','N','U','E','D','S','E'};
  /**
    // Initializes a random 4-by-4 Boggle board.
    // (by rolling the Hasbro dice)
   */
    public BoggleBoard() {
      board = new char[4][4];
      int[] indices = new int[rows()*cols()];
      for (int i = 0; i < indices.length; i++) {
        indices[i] = i;
      }
      StdRandom.shuffle(indices);
      for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < cols(); j++) {
          board[i][j] = hasbroDice[indices[i*rows() + j]];
        }
      }
    }

    /**
    // Initializes a random m-by-n Boggle board.
    // (using the frequency of letters in the English language)
     */
    public BoggleBoard(int m, int n) {
      board = new char[m][n];
      throw new NotImplementedException();
    }

    /**
    // Initializes a Boggle board from the specified filename.
     */
    public BoggleBoard(String filename) {
      In in = new In(filename);
      int nrow = in.readInt();
      int ncol = in.readInt();
      String[] matrix = in.readAllStrings();
      board = new char[nrow][ncol];
      for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
          board[i][j] = matrix[i*rows() + j].charAt(0);
        }
      }

    }

    /**
    // Initializes a Boggle board from the 2d char array.
    // (with 'Q' representing the two-letter sequence "Qu")
     */
    public BoggleBoard(char[][] a) {
      board = new char[a.length][a[0].length];
      for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < cols(); j++) {
          board[i][j] = a[i][j];
        }
      }
    }

    /**
    // Returns the number of rows.
     */
    public int rows() {
      return board.length;
    }

    /**
    // Returns the number of columns.
     */
    public int cols() {
      return board[0].length;
    }

    /**
    // Returns the letter in row i and column j.
    // (with 'Q' representing the two-letter sequence "Qu")
     */
    public char getLetter(int i, int j) {
      return board[i][j];
    }

    /**
    // Returns a string representation of the board.
     */
    public String toString() {
      StringBuilder sb = new StringBuilder();
      for (int i = 0; i < board.length; i++) {
        for (int j = 0; j < board.length; j++) {
          sb.append(board[i][j]);
          sb.append(" ");
        }
        sb.append("\n");
      }
      return sb.toString();
    }
}
