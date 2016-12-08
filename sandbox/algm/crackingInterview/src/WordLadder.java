import java.util.*;

public class WordLadder {
        /**
          * @param start, a string
          * @param end, a string
          * @param dict, a set of string
          * @return an integer
          */
    //optimize for large dictionary
    public int ladderLength(String start, String end, Set<String> dict) {
                // write your code here
            if (start.equals(end)) {
                return 2;
    } else if (dict.size() == 0) {
        return 0;
    }
    //bfs
    Deque<String> toVisit = new LinkedList<String> ();
    int lvl = 0;
    Set<String> dictCopy = cloneDict(dict); //keep original dict
    //Set has no clone method
    toVisit.offerLast(start);
    assert(!dictCopy.contains(end)); //assume end is not in dict
    dictCopy.add(end);
    dictCopy.remove(start);
    while (!toVisit.isEmpty()) {
        lvl++;
        int currentLevelSize = toVisit.size();
        for (int i = 0; i < currentLevelSize; i++) {
            String cur = toVisit.pollFirst();
            if (cur.equals(end)) {          
                return lvl;
            }
            addChildren(cur, toVisit, dictCopy);
        }
    }
    return 0;
            }
        private Set<String> cloneDict(Set<String> dict) {
            Set<String> copy = new HashSet<String> ();
            for (String i : dict) {
                copy.add(i);
            }
            return copy;
        }
        /*check and add all unvisited children that are in dictCopy to toVisit */
    private void addChildren(String cur, Deque<String> toVisit,  Set<String> dictCopy) {
            for (int i = 0; i < cur.length(); i++) {
                for (char c = 'a'; c <= 'z'; c++) {
                    String newStr = cur.substring(0, i) + c + cur.substring(i + 1, cur.length());
                    if (dictCopy.contains(newStr)) {
                        toVisit.offerLast(newStr);
                        dictCopy.remove(newStr);
                    }
                }
            }
    }
    public static void main(String[] args) {
        Set<String> dict = new HashSet<String> ();
        dict.add("hot");
        dict.add("dot");
        dict.add("dog");
        dict.add("lot");
        dict.add("log");
        WordLadder test = new WordLadder();
        System.out.println(test.ladderLength("hit", "cog", dict));
    }
}
