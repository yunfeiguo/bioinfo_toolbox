import java.util.*;
public class WordLadderII {
    //BFS, iterative
    public List<List<String>> findLadders(String start, String end, Set<String> dict){
        if(dict == null){
            return new ArrayList<List<String>>();
        }
        dict.add(end);
        // bfs
        List<List<String>> result = new ArrayList<List<String>>();
        Deque<List<String>> queue = new LinkedList<List<String>>();
        List<String> startingState = new ArrayList<String>();
        Set<String> visited = new HashSet<String>();
        visited.add(start);
        startingState.add(start);
        queue.offerFirst(startingState);
        while(!queue.isEmpty()){
            int levelSize = queue.size();
            boolean flag = false; // test whether or not we have reached end
            Set<String> set = new HashSet<>();
            for(int i = 0; i < levelSize; i++){
                List<String> current = queue.pollLast();
                List<List<String>> nextStrings = getNextStrings(current, dict);
                for(List<String> next : nextStrings){
                    if(visited.contains(next.get(next.size() - 1))){
                        continue;
                    }
                    if(!flag && next.get(next.size() - 1).equals(end)){
                        flag = true;
                        result.add(next);
                    }else if(flag && next.get(next.size() - 1).equals(end)){
                        result.add(next);
                    }else{
                        queue.offerFirst(next);
                    }
                    set.add(next.get(next.size() - 1));
                }
            }
            for(String string : set){
                visited.add(string);
            }
            if(flag){
                return result;
            }
        }
        return result;
    }
    private List<List<String>> getNextStrings(List<String> current, Set<String> dict){
        List<List<String>> next = new ArrayList<List<String>>();
        String currentLast = current.get(current.size() - 1);
        for(String dictString : dict){
            if(differOne(dictString, currentLast)){
                current.add(dictString);
                next.add(new ArrayList(current));
                current.remove(current.size() - 1);
            }
        }
        return next;
    }
    private boolean differOne(String dictString, String current){
        int index = 0;
        for(int i = 0; i < dictString.length(); i++){
            if(dictString.charAt(i) != current.charAt(i)){
                index++;
            }
        }
        return index == 1;
    }   
    public static void main(String[] args) {
        WordLadderII test = new WordLadderII();
        Set<String> dict = new HashSet<String>();
        //String[] allString = {"si","go","se","cm","so","ph","mt","db","mb","sb","kr","ln","tm","le","av","sm","ar","ci","ca","br","ti","ba","to","ra","fa","yo","ow","sn","ya","cr","po","fe","ho","ma","re","or","rn","au","ur","rh","sr","tc","lt","lo","as","fr","nb","yb","if","pb","ge","th","pm","rb","sh","co","ga","li","ha","hz","no","bi","di","hi","qa","pi","os","uh","wm","an","me","mo","na","la","st","er","sc","ne","mn","mi","am","ex","pt","io","be","fm","ta","tb","ni","mr","pa","he","lr","sq","ye"};
        String[] allString = {"hot", "cog", "dog", "tot", "hog", "hop", "pot","dot"};
        for (String s : allString) {
            dict.add(s);
        }
        test.findLadders("hot", "dog", dict);
    }

}
