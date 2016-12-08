/*
 * String is literal, immutable
 */
public class StringPractice {
    public static void main(String[] args) {
        char[] abc = {'a','b', 'c'};
        String str = new String(abc);
        System.out.println(str);
        String str2 = "abc";
        System.out.println(str == str2);
        System.out.println(str.equals(str2));
        str.concat(str2);
        System.out.println(str.concat(str2));
        String[] split = str.split("");
        StringBuilder join = new StringBuilder();
        for (String s : split) {
            join.append(s);
            join.append("*");
        }
        join.delete(join.length() - 1, join.length());
        System.out.println(join.toString());
    }   

}
