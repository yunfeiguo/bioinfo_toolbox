package yunfeiImplementAlgs4;

import edu.princeton.cs.algs4.StdIn;
import edu.princeton.cs.algs4.StdOut;

/* output the longest repeated
 * string from stdin
 */
public class LRS {
	public static void main(String[] args) {
		String text = StdIn.readAll();
		//StdOut.println("input:"+text);
		int N = text.length();
		SuffixArray sa = new SuffixArray(text);
		String lrs = "";
		for (int i = 1; i < N; i++){
			int length = sa.lcp(i);
			//StdOut.println(Integer.toString(i));
			//StdOut.println(sa.select(i));
			if (length > lrs.length())
				lrs = sa.select(i).substring(0, length);			
		}
		StdOut.println(lrs);
	}
}
