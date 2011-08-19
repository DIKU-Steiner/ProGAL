package ProGAL.proteins.structure;

public class CBond {
	public final Atom a1, a2;
	
	public CBond(Atom a1, Atom a2){
		this.a1=a1;
		this.a2=a2;
	}

	public double length() {
		return a1.distance(a2);
	}
	
}
