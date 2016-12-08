
public class ConvertDecimal5_2 {
	public static void main(String[] args) {
		String s = "1.5";
		//String s2 = ".5";
		//System.out.print(Double.parseDouble(s2));
		System.out.println(ConvertDecimal5_2.convert(s).equals("1.1"));
	}
	public static String convert(String s) {
		int maxFrac = 32; //max length of fractional part
		int dot = s.indexOf('.');		
		if (dot == -1) {
			//no fractional part
			return(Integer.toBinaryString(Integer.parseInt(s)));		
		}
		int intPart = Integer.parseInt(s.substring(0, dot));
		double fracPart = Double.parseDouble(s.substring(dot));
		String fracRet = "";
		for (int i = 0; i < maxFrac; i++) {
			fracPart = fracPart * 2;
			short cur = (short) Math.floor(fracPart);
			fracPart = fracPart - cur;
			fracRet = fracRet + cur;
			if(fracPart == 0) break;						
		}
		if (fracPart != 0) {
			return("ERROR");
		}		
		return(Integer.toBinaryString(intPart) + "." + fracRet);
	}

}
