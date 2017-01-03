package yunfeiImplementAlgs4;

import java.io.File;

import edu.princeton.cs.algs4.Insertion;

public class ListFiles {
    public static void main(String[] args) {
        File dir = new File(args[0]);
        File[] files = dir.listFiles();
        Insertion.sort(files);
        for (File s : files)
            System.out.println(s.getName());
    }

}
