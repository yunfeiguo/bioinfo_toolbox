import java.util.HashMap;
import java.util.Map;
public class StrStr {
	    /**
	     * Returns a index to the first occurrence of target in source,
	     * or -1  if target is not part of source.
	     * @param source string to be scanned.
	     * @param target string containing the sequence of characters to match.
	     */
	    public static int strStr(String source, String target) {
	        //write your code here
	        //from 0 in source, check if target is a substring in source at 0
	        //if length from 0 is shorter than target, return -1 directly
	        int sLen = source.length();
	        int tLen = target.length();
	        if(tLen == 0) return(0);//special case

	        for(int i = 0; i < sLen; i++) {
	            if (sLen - i < tLen ) {
	                return(-1);//remaining length shorter than target
	            }
	            System.out.println(source.substring(i, i + tLen));
	            if (source.substring(i, i + tLen).equals(target)) {
	                
	                return(i);
	            }
	        }
	        return -1;
	    }
public static void main(String[] args) {
	StrStr.strStr("abcdabcdefg", "bcd");
	HashMap<String, Integer> h = new HashMap<String, Integer>();
}
}
