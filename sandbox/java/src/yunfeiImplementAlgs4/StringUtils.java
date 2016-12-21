package yunfeiImplementAlgs4;

public class StringUtils {
	public static String join(String[] s, String sep) {
		String ret = new String();
		for (int i = 0; i < s.length; i++) {			
			if (i > 0) {
				ret = sep + ret;
			}
			ret = ret + s[i];
		}
		return(ret);
	}

}
