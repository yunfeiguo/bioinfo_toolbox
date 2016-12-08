public class Solution {
    /**
     * @param strs: A list of strings
     * @return: The longest common prefix
     */
    public String longestCommonPrefix(String[] strs) {
        // write your code here
        if (strs == null || strs.length == 0) return(new String());
        int l = strs[0].length();
        for (String i : strs) {
            l = Math.min(l, i.length());
        }
        //iterate over 0 to end of shortest string
        //stop any at position where there are different char
        //across all strings
        int i = 0;
        boolean stop = false;
        while (i < l) {
            char c = strs[0].charAt(i);
            for (String s: strs) {
                if (s.charAt(i) != c) {
                    stop = true;
                    break;
                }
            }
            if (stop) break;
            i++;
        }
        return(strs[0].substring(0, i));
    }
}
