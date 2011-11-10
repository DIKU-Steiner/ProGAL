package ProGAL.proteins.belta;

public class PrimaryStructure {
	/** The sequence of a protein chain */
	public final String sequence;

	public PrimaryStructure(String seq){
		this.sequence = seq;
	}
	
	public String toString(){
		return sequence;
	}
}
