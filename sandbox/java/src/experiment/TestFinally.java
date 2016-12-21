package experiment;

import java.io.File;

/**
 * Created by guoy28 on 10/4/16.
 */
public class TestFinally {
    public static void main(String[] args) {
        try {
            File f = new File("/tmp/1.txt");
            return;
        } catch (Exception e) {
            System.out.println(e.getMessage());
        } finally {
            System.out.println("finally");
        }
    }
}
