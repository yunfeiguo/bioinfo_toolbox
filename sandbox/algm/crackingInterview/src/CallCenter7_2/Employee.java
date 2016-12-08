package CallCenter7_2;

public class Employee {
	private boolean available;
	//will be used to determine whether
	//this person can handle a call
	private int rank;	
	private int id;
	private String name;
	public Employee(int r, int i, String n) {
		rank = r;
		setId(i);
		setName(n);
		setAvailable(true);
	}
	public boolean isAvailable() {
		return(available);
	}	
	public int getRank() {
		return(rank);
	}
	public void setRank(int r) {
		rank = r;
	}
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}
	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}
	/**
	 * @param id the id to set
	 */
	public void setId(int id) {
		this.id = id;
	}
	public boolean canHandleCall(PhoneCall call) {
		System.out.println(rank + " " + call.getRank());
		return(rank >= call.getRank());
	}
	/**if can handle call, set available to false, return true
	 * otherwise, return false
	 */
	public boolean handleCall(PhoneCall call) {
		if(isAvailable() && canHandleCall(call)) {
			setAvailable(false);
			return(true);
		} else {
			return(false);
		}
	}
	public void setAvailable(boolean b) {
		available = b;
	}
}
