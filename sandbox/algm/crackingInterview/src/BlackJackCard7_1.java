
public class BlackJackCard7_1 extends Card7_1{
	public BlackJackCard7_1(int r, Suit s) {super(r, s);}
	public int value() {
		int r = super.value();
		if (r == 1) return(11);
		if (r < 10) return(r);
		return(10);
	}
	boolean isAce() {
		return(super.value() == 1);
	}
}
