package CallCenter7_2;

public class CallCenter {
	private Fresher[] freshers;//we could use ArrayList for all employees
	//because employees may come and go
	private TechnicalLead TL;
	private ProjectManager PM;
	public CallCenter(int n) {
		freshers = new Fresher[n];
		for (int i = 0; i < freshers.length; i++) {
			freshers[i] = new Fresher(i, "Fresher"+i);
		}
		TL = new TechnicalLead(n, "TL");
		PM = new ProjectManager(n+1, "PM");
	}
	public static void main(String[] args) {
		//initiate the operation
		CallCenter cc = new CallCenter(10);
		PhoneCall call = new PhoneCall(1);//we can use a queue for all the incoming calls
		Employee e = cc.getCallHandler(call);
		System.out.println(e == null? "Nobody": e.getName());
		PhoneCall call2 = new PhoneCall(1);
		Employee e2 = cc.getCallHandler(call2);
		System.out.println(e2 == null? "Nobody": e2.getName());
	}
	//we can have a class for call handler
	public Employee getCallHandler(PhoneCall call) {
		//look for available employees
		//assign proper employee to this call
		for(Fresher i : freshers) {
			if(i.handleCall(call)) {
				return(i);
			}
		}
		if(TL.handleCall(call)) {
			return(TL);
		}
		if(PM.handleCall(call)) {
			return(PM);
		}
		return(null);
	}

}
